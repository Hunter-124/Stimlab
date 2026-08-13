#include "modules/docking/VinaParser.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace biocad::docking {
namespace {

namespace chem = biocad::chem;

bool startsWith(const std::string& s, const char* prefix) {
    return s.rfind(prefix, 0) == 0;
}

// Parse the three doubles after "REMARK VINA RESULT:". Tolerates extra whitespace
// and missing rmsd fields (defaults 0). Returns false if no affinity is present.
bool parseVinaResult(const std::string& line, double& aff, double& lb, double& ub) {
    const std::string tag = "VINA RESULT:";
    const auto pos = line.find(tag);
    if (pos == std::string::npos) return false;
    std::istringstream is(line.substr(pos + tag.size()));
    aff = lb = ub = 0.0;
    if (!(is >> aff)) return false;
    is >> lb;  // optional
    is >> ub;  // optional
    return std::isfinite(aff);
}

// Extract x,y,z from a PDB(QT) ATOM/HETATM record using fixed columns
// (31-38, 39-46, 47-54, 1-based). Falls back to whitespace tokenisation if the
// fixed-column slice is not numeric (some writers shift columns slightly).
bool parseAtomXyz(const std::string& line, double& x, double& y, double& z) {
    auto slice = [&](size_t start, size_t len, double& out) -> bool {
        if (line.size() < start + len) return false;
        try {
            out = std::stod(line.substr(start, len));
            return std::isfinite(out);
        } catch (...) {
            return false;
        }
    };
    if (slice(30, 8, x) && slice(38, 8, y) && slice(46, 8, z)) return true;

    // Fallback: tokenise and grab the first three finite doubles after column 6.
    std::istringstream is(line.size() > 6 ? line.substr(6) : std::string());
    std::vector<double> nums;
    std::string tok;
    while (is >> tok) {
        try {
            double v = std::stod(tok);
            if (std::isfinite(v)) nums.push_back(v);
        } catch (...) {
        }
        if (nums.size() >= 5) break;  // x,y,z are early; don't scan the whole line
    }
    if (nums.size() >= 3) { x = nums[0]; y = nums[1]; z = nums[2]; return true; }
    return false;
}

// Build a pose conformer: parsed coordinates wear the reference molecule's element
// list, bonds and heavyCount. If the parsed atom count differs from the reference
// (e.g. a writer dropped/added an H), we keep the reference z/bonds for the atoms
// we have and clamp heavyCount so indices stay valid for the 3D viewer.
chem::Conformer makeConformer(const std::vector<chem::Vec3>& coords,
                              const chem::Conformer& reference,
                              const std::vector<int>* order) {
    // Order-aware path: the writer emitted atoms in a permuted (torsion-tree) order,
    // so scatter each output coordinate back onto its source atom in `reference`. The
    // full topology (z, bonds, heavyCount) is preserved; atoms not emitted (nonpolar
    // H) keep their embedded positions.
    if (order && !order->empty()) {
        chem::Conformer c = reference;
        const int n = static_cast<int>(coords.size());
        for (int k = 0; k < n && k < static_cast<int>(order->size()); ++k) {
            const int idx = (*order)[k];
            if (idx >= 0 && idx < static_cast<int>(c.pos.size())) c.pos[idx] = coords[k];
        }
        return c;
    }

    chem::Conformer c;
    c.pos = coords;
    const int n = static_cast<int>(coords.size());
    c.z.resize(n, 6);
    for (int i = 0; i < n && i < reference.size(); ++i) c.z[i] = reference.z[i];
    // Keep only bonds whose endpoints exist in this pose.
    for (const auto& b : reference.bonds)
        if (b.first < n && b.second < n) c.bonds.push_back(b);
    c.heavyCount = std::min(reference.heavyCount, n);
    return c;
}

}  // namespace

std::vector<DockPose> parseVinaPdbqt(const std::string& output,
                                     const chem::Conformer& reference,
                                     const std::vector<int>* serialToConf) {
    std::vector<DockPose> poses;
    std::istringstream in(output);
    std::string line;

    bool inModel = false;
    bool haveResult = false;
    double aff = 0.0, lb = 0.0, ub = 0.0;
    std::vector<chem::Vec3> coords;

    auto flush = [&]() {
        // Emit a pose only if we saw a VINA RESULT and at least one atom.
        if (haveResult && !coords.empty()) {
            DockPose p;
            p.affinityKcalPerMol = aff;
            p.rmsdLb = lb;
            p.rmsdUb = ub;
            p.ligand = makeConformer(coords, reference, serialToConf);
            poses.push_back(std::move(p));
        }
        inModel = false;
        haveResult = false;
        aff = lb = ub = 0.0;
        coords.clear();
    };

    while (std::getline(in, line)) {
        // Normalise a trailing CR from Windows-authored output.
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (startsWith(line, "MODEL")) {
            if (inModel) flush();  // tolerate a missing ENDMDL between MODELs
            inModel = true;
            haveResult = false;
            coords.clear();
            continue;
        }
        if (startsWith(line, "ENDMDL")) {
            flush();
            continue;
        }
        if (line.find("VINA RESULT:") != std::string::npos) {
            double a, l, u;
            if (parseVinaResult(line, a, l, u)) {
                aff = a; lb = l; ub = u; haveResult = true;
            }
            continue;
        }
        if (startsWith(line, "ATOM") || startsWith(line, "HETATM")) {
            double x, y, z;
            if (parseAtomXyz(line, x, y, z)) coords.push_back({x, y, z});
            continue;
        }
    }
    // A final MODEL without a closing ENDMDL still counts.
    if (inModel) flush();

    // Rank most-negative (strongest) first and renumber.
    std::sort(poses.begin(), poses.end(), [](const DockPose& a, const DockPose& b) {
        return a.affinityKcalPerMol < b.affinityKcalPerMol;
    });
    for (int i = 0; i < static_cast<int>(poses.size()); ++i) poses[i].rank = i + 1;
    return poses;
}

}  // namespace biocad::docking
