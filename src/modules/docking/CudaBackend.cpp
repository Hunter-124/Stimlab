// modules/docking/CudaBackend.cpp - CUDA GPU docking backend (host side).
//
// Plain C++ (MSVC-compiled): parses the prepared rigid receptor PDBQT, generates a grid
// of rigid ligand poses (rotations x translations within the receptor box), hands them
// to the GPU scorer (CudaScore.cu), and returns the best-scoring poses. Only the small
// POD scoring kernel is compiled by nvcc; all the chemistry/IO lives here so nvcc never
// sees the C++20 headers. Deterministic (a fixed pose grid - no RNG), mirroring Vina's
// reproducible --seed 1.
#include "modules/docking/CudaBackend.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "modules/docking/CudaScore.h"

namespace stimlab::docking {
namespace {

constexpr double kPi = 3.14159265358979323846;

// Vina xs surface radii (Angstrom) + the 2-class type (0 = hydrophobic, 1 = polar).
void ligTypeRadius(int z, int& type, float& rad) {
    switch (z) {
        case 7:  type = 1; rad = 1.8f; break;  // N
        case 8:  type = 1; rad = 1.7f; break;  // O
        case 16: type = 0; rad = 2.0f; break;  // S
        case 15: type = 0; rad = 2.1f; break;  // P
        case 9:  type = 0; rad = 1.5f; break;  // F
        case 17: type = 0; rad = 1.8f; break;  // Cl
        case 35: type = 0; rad = 2.0f; break;  // Br
        case 53: type = 0; rad = 2.2f; break;  // I
        default: type = 0; rad = 1.9f; break;  // C and everything else: carbon-like
    }
}

// Map an AutoDock4 receptor atom type token (C, A, N, NA, OA, SA, S, P, ...) to the
// same 2-class type + radius used for the ligand.
void recTypeRadius(const std::string& t, int& type, float& rad) {
    if (t == "O" || t == "OA") { type = 1; rad = 1.7f; return; }
    if (t == "N" || t == "NA") { type = 1; rad = 1.8f; return; }
    if (t == "S" || t == "SA") { type = 0; rad = 2.0f; return; }
    if (t == "P")              { type = 0; rad = 2.1f; return; }
    if (t == "Cl" || t == "CL"){ type = 0; rad = 1.8f; return; }
    if (t == "Br" || t == "BR"){ type = 0; rad = 2.0f; return; }
    if (t == "F")              { type = 0; rad = 1.5f; return; }
    if (t == "I")              { type = 0; rad = 2.2f; return; }
    type = 0; rad = 1.9f;  // C / A (aromatic) / unknown -> hydrophobic carbon-like
}

struct RecAtoms {
    std::vector<float> xyz;   // [x0,y0,z0, ...]
    std::vector<int>   type;
    std::vector<float> rad;
    int count = 0;
};

// Parse a prepared rigid receptor PDBQT, keeping only atoms within (box half-extent +
// 8 A cutoff) of the box center so the GPU work is bounded to the pocket neighborhood.
RecAtoms parseReceptor(const std::string& path, const DockBox& box) {
    RecAtoms r;
    std::ifstream in(path, std::ios::binary);
    if (!in) return r;
    const double keep = std::max({box.sx, box.sy, box.sz}) * 0.5 + 8.0;
    const double keep2 = keep * keep;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 54) continue;
        const std::string rec = line.substr(0, 6);
        if (rec != "ATOM  " && rec != "HETATM") continue;
        double x, y, z;
        try {
            x = std::stod(line.substr(30, 8));
            y = std::stod(line.substr(38, 8));
            z = std::stod(line.substr(46, 8));
        } catch (...) {
            continue;
        }
        const double dx = x - box.cx, dy = y - box.cy, dz = z - box.cz;
        if (dx * dx + dy * dy + dz * dz > keep2) continue;  // outside the pocket window
        // AutoDock type = last whitespace-separated token on the line.
        std::istringstream ss(line);
        std::string tok, last;
        while (ss >> tok) last = tok;
        int t;
        float rad;
        recTypeRadius(last, t, rad);
        r.xyz.push_back(static_cast<float>(x));
        r.xyz.push_back(static_cast<float>(y));
        r.xyz.push_back(static_cast<float>(z));
        r.type.push_back(t);
        r.rad.push_back(rad);
    }
    r.count = static_cast<int>(r.type.size());
    return r;
}

// Row-major 3x3 rotation R = Rz(a) Ry(b) Rx(g).
void eulerR(double a, double b, double g, double R[9]) {
    const double ca = std::cos(a), sa = std::sin(a);
    const double cb = std::cos(b), sb = std::sin(b);
    const double cg = std::cos(g), sg = std::sin(g);
    R[0] = ca * cb; R[1] = ca * sb * sg - sa * cg; R[2] = ca * sb * cg + sa * sg;
    R[3] = sa * cb; R[4] = sa * sb * sg + ca * cg; R[5] = sa * sb * cg - ca * sg;
    R[6] = -sb;     R[7] = cb * sg;                R[8] = cb * cg;
}

}  // namespace

bool CudaBackend::available() const { return cudaDockAvailable(); }

DockJobResult CudaBackend::dock(const chem::Molecule& graph, const chem::Conformer& ligand3d,
                                const ReceptorTarget& target) const {
    DockJobResult res;
    res.engine = "CUDA GPU (rigid, Vina scoring)";
    res.real = false;
    res.targetId = target.id;

    if (!cudaDockAvailable()) { res.log = "No CUDA device available."; return res; }
    if (ligand3d.empty() || ligand3d.heavyCount <= 0) {
        res.log = "Ligand has no 3D conformer.";
        return res;
    }
    if (target.receptorPath.empty()) {
        res.log = "No prepared receptor for this target (provision it to enable GPU docking).";
        return res;
    }
    const RecAtoms rec = parseReceptor(target.receptorPath, target.box);
    if (rec.count == 0) {
        res.log = "Prepared receptor had no atoms inside the box window.";
        return res;
    }

    // Ligand heavy atoms (original conformer coordinates) + types/radii + centroid.
    const int nLig = ligand3d.heavyCount;
    std::vector<float> ligXYZ(static_cast<size_t>(nLig) * 3);
    std::vector<int>   ligType(nLig);
    std::vector<float> ligRad(nLig);
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < nLig; ++i) {
        const auto& p = ligand3d.pos[i];
        ligXYZ[i * 3 + 0] = static_cast<float>(p.x);
        ligXYZ[i * 3 + 1] = static_cast<float>(p.y);
        ligXYZ[i * 3 + 2] = static_cast<float>(p.z);
        const int z = (i < static_cast<int>(ligand3d.z.size())) ? ligand3d.z[i] : 6;
        ligTypeRadius(z, ligType[i], ligRad[i]);
        cx += p.x; cy += p.y; cz += p.z;
    }
    cx /= nLig; cy /= nLig; cz /= nLig;

    // Pose grid: rotations x a translation grid spanning the box (<= ~12 steps/axis so
    // the pose count stays bounded for any box size). Each pose is a 3x4 transform that
    // rotates the ligand about its centroid then places that centroid at a box point.
    std::vector<std::array<double, 9>> rots;
    const int kYaw = 8, kPitch = 4, kRoll = 2;  // 64 orientations
    for (int a = 0; a < kYaw; ++a)
        for (int b = 0; b < kPitch; ++b)
            for (int g = 0; g < kRoll; ++g) {
                std::array<double, 9> R{};
                eulerR(2.0 * kPi * a / kYaw, kPi * b / kPitch, kPi * g / kRoll, R.data());
                rots.push_back(R);
            }

    auto axisGrid = [](double center, double size) {
        const double step = std::max(2.0, size / 11.0);  // ~<=12 positions per axis
        std::vector<double> v;
        for (double t = center - size * 0.5 + step * 0.5; t <= center + size * 0.5; t += step)
            v.push_back(t);
        if (v.empty()) v.push_back(center);
        return v;
    };
    const std::vector<double> txs = axisGrid(target.box.cx, target.box.sx);
    const std::vector<double> tys = axisGrid(target.box.cy, target.box.sy);
    const std::vector<double> tzs = axisGrid(target.box.cz, target.box.sz);

    const size_t nPoses = rots.size() * txs.size() * tys.size() * tzs.size();
    std::vector<float> xform(nPoses * 12);
    size_t pi = 0;
    for (const auto& R : rots) {
        // R * centroid (constant per rotation).
        const double rcx = R[0] * cx + R[1] * cy + R[2] * cz;
        const double rcy = R[3] * cx + R[4] * cy + R[5] * cz;
        const double rcz = R[6] * cx + R[7] * cy + R[8] * cz;
        for (double tx : txs)
            for (double ty : tys)
                for (double tz : tzs) {
                    float* M = &xform[pi * 12];
                    M[0] = (float)R[0]; M[1] = (float)R[1]; M[2] = (float)R[2]; M[3] = (float)(tx - rcx);
                    M[4] = (float)R[3]; M[5] = (float)R[4]; M[6] = (float)R[5]; M[7] = (float)(ty - rcy);
                    M[8] = (float)R[6]; M[9] = (float)R[7]; M[10] = (float)R[8]; M[11] = (float)(tz - rcz);
                    ++pi;
                }
    }

    std::vector<float> scores(nPoses);
    if (!cudaScorePoses(rec.xyz.data(), rec.type.data(), rec.rad.data(), rec.count, ligXYZ.data(),
                        ligType.data(), ligRad.data(), nLig, xform.data(),
                        static_cast<int>(nPoses), scores.data())) {
        res.log = "GPU scoring failed.";
        return res;
    }

    // Rank: take the most-favorable (most negative) poses. partial_sort the indices.
    const int kKeep = std::min<int>(9, static_cast<int>(nPoses));
    std::vector<int> idx(nPoses);
    std::iota(idx.begin(), idx.end(), 0);
    std::partial_sort(idx.begin(), idx.begin() + kKeep, idx.end(),
                      [&](int a, int b) { return scores[a] < scores[b]; });

    for (int k = 0; k < kKeep; ++k) {
        const float* M = &xform[static_cast<size_t>(idx[k]) * 12];
        // Build the docked conformer: apply this pose's transform to EVERY atom (heavy
        // + H) so the 3D viewer shows the full ligand in the pocket.
        chem::Conformer out;
        out.z = ligand3d.z;
        out.bonds = ligand3d.bonds;
        out.heavyCount = ligand3d.heavyCount;
        out.pos.resize(ligand3d.pos.size());
        for (size_t a = 0; a < ligand3d.pos.size(); ++a) {
            const auto& p = ligand3d.pos[a];
            out.pos[a].x = M[0] * p.x + M[1] * p.y + M[2] * p.z + M[3];
            out.pos[a].y = M[4] * p.x + M[5] * p.y + M[6] * p.z + M[7];
            out.pos[a].z = M[8] * p.x + M[9] * p.y + M[10] * p.z + M[11];
        }
        DockPose pose;
        pose.rank = k + 1;
        pose.affinityKcalPerMol = scores[idx[k]];
        pose.rmsdLb = 0.0;
        pose.rmsdUb = 0.0;
        pose.ligand = std::move(out);
        res.poses.push_back(std::move(pose));
    }

    res.real = true;
    res.log = "Docked on GPU: " + std::to_string(nPoses) + " rigid poses scored (Vina function), " +
              std::to_string(rec.count) + " receptor atoms.";
    (void)graph;
    return res;
}

}  // namespace stimlab::docking
