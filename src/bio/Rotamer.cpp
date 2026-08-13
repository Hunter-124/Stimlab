#include "bio/Rotamer.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>

#include "bio/Score.h"
#include "core/Error.h"
#include "packs/Pack.h"

namespace biocad::bio {
namespace {

struct Vec3 {
    double x = 0, y = 0, z = 0;
};

Vec3 sub(const Point3d& a, const Point3d& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double len(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 scale(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
Vec3 unit(const Vec3& a) {
    const double n = len(a);
    return n > 0 ? scale(a, 1.0 / n) : a;
}

constexpr double kDeg = 3.14159265358979323846 / 180.0;

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t");
    return s.substr(b, e - b + 1);
}

}  // namespace

double dihedralDegrees(const Point3d& p0, const Point3d& p1, const Point3d& p2,
                       const Point3d& p3) {
    const Vec3 b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
    const Vec3 u = unit(b1);
    const Vec3 v{b0.x - dot(b0, u) * u.x, b0.y - dot(b0, u) * u.y, b0.z - dot(b0, u) * u.z};
    const Vec3 w{b2.x - dot(b2, u) * u.x, b2.y - dot(b2, u) * u.y, b2.z - dot(b2, u) * u.z};
    return std::atan2(dot(cross(u, v), w), dot(v, w)) / kDeg;
}

Point3d placeAtom(const Point3d& a, const Point3d& b, const Point3d& c, double bond,
                  double angleDeg, double torsionDeg) {
    // Natural extension reference frame. The local frame is built on c: bc is the
    // chain direction and n the normal of the a-b-c plane, taken as (b - a) x bc.
    // The orientation of n is not cosmetic: taking (a - b) x bc instead flips both
    // in-plane axes and silently places every atom at torsion + 180 degrees, which
    // the placeAtom/dihedralDegrees round-trip test catches.
    const Vec3 bc = unit(sub(c, b));
    const Vec3 n = unit(cross(sub(b, a), bc));
    const Vec3 m = cross(n, bc);

    const double ang = angleDeg * kDeg, tor = torsionDeg * kDeg;
    // d in the local frame: the bond makes (180 - angle) with bc.
    const double d0 = -bond * std::cos(ang);
    const double d1 = bond * std::sin(ang) * std::cos(tor);
    const double d2 = bond * std::sin(ang) * std::sin(tor);
    return {c.x + d0 * bc.x + d1 * m.x + d2 * n.x,
            c.y + d0 * bc.y + d1 * m.y + d2 * n.y,
            c.z + d0 * bc.z + d1 * m.z + d2 * n.z};
}

// ---------------------------------------------------------------------------
// Library
// ---------------------------------------------------------------------------

const ResidueRotamers* RotamerLibrary::residue(const std::string& threeLetter) const {
    for (const auto& r : residues_)
        if (r.name == threeLetter) return &r;
    return nullptr;
}

std::vector<RotamerEntry> RotamerLibrary::rotamersAt(const std::string& threeLetter, double phi,
                                                     double psi,
                                                     bool& usedBackboneIndependent) const {
    usedBackboneIndependent = true;
    const ResidueRotamers* r = residue(threeLetter);
    if (!r) return {};
    const double w = static_cast<double>(binDegrees_);
    const int pb = static_cast<int>(std::floor((phi + 180.0) / w)) * binDegrees_ - 180;
    const int sb = static_cast<int>(std::floor((psi + 180.0) / w)) * binDegrees_ - 180;
    for (const auto& bin : r->bins) {
        if (bin.phi == pb && bin.psi == sb) {
            usedBackboneIndependent = false;
            return bin.rotamers;
        }
    }
    return r->backboneIndependent;
}

RotamerLibrary parseRotamerLibrary(const nlohmann::json& j) {
    if (!j.is_object()) throw Error::parse("rotamer pack: the document is not an object");
    const int version = j.value("schemaVersion", 0);
    if (version != 1)
        throw Error::parse("rotamer pack: schemaVersion " + std::to_string(version) +
                           " is not supported by this build (expected 1)");

    RotamerLibrary lib;
    lib.name_ = j.value("libraryName", std::string{});
    lib.method_ = j.value("method", std::string{});
    if (lib.name_.empty())
        throw Error::parse("rotamer pack: libraryName is required - a rebuild with no named "
                           "library is not reproducible");
    if (j.contains("attribution"))
        for (const auto& a : j.at("attribution")) lib.attribution_.push_back(a.get<std::string>());
    if (j.contains("notDunbrack")) lib.attribution_.push_back(j.at("notDunbrack").get<std::string>());
    lib.binDegrees_ = j.value("phiPsiBinDegrees", 60);
    if (j.contains("dataset")) {
        const auto& d = j.at("dataset");
        lib.datasetEntries_ = d.value("entryCount", 0);
        lib.residuesMeasured_ = d.value("residuesMeasured", 0);
    }

    const auto readRotamers = [](const nlohmann::json& arr) {
        std::vector<RotamerEntry> out;
        for (const auto& e : arr) {
            RotamerEntry r;
            r.name = e.value("name", std::string{});
            r.count = e.value("count", 0);
            r.probability = e.value("probability", 0.0);
            if (e.contains("chi")) r.chi = e.at("chi").get<std::vector<double>>();
            if (e.contains("chiSd")) r.chiSd = e.at("chiSd").get<std::vector<double>>();
            out.push_back(std::move(r));
        }
        return out;
    };

    if (!j.contains("residues")) throw Error::parse("rotamer pack: no `residues` object");
    for (const auto& [name, r] : j.at("residues").items()) {
        ResidueRotamers rr;
        rr.name = name;
        rr.chiCount = r.value("chiCount", 0);
        rr.nonRotamericChi = r.value("nonRotamericChi", 0);
        rr.symmetricChi = r.value("symmetricChi", 0);
        rr.backboneIndependentCount = r.value("backboneIndependentCount", 0);
        for (const auto& b : r.value("build", nlohmann::json::array())) {
            BuildAtom a;
            a.atom = b.at("atom").get<std::string>();
            const auto p = b.at("parents");
            if (p.size() != 3)
                throw Error::parse("rotamer pack: " + name + "/" + a.atom +
                                   " needs exactly three parents");
            for (int i = 0; i < 3; ++i) a.parents[i] = p[static_cast<std::size_t>(i)];
            a.bondLength = b.value("bondLength", 0.0);
            a.bondAngle = b.value("bondAngle", 0.0);
            a.torsion = b.value("torsion", 0.0);
            a.torsionSd = b.value("torsionSd", 0.0);
            a.chi = b.value("chi", 0);
            a.observations = b.value("observations", 0);
            if (a.bondLength <= 0.0)
                throw Error::parse("rotamer pack: " + name + "/" + a.atom +
                                   " has no measured bond length");
            rr.build.push_back(std::move(a));
        }
        if (r.contains("backboneIndependent"))
            rr.backboneIndependent = readRotamers(r.at("backboneIndependent"));
        for (const auto& b : r.value("bins", nlohmann::json::array())) {
            PhiPsiBin bin;
            bin.phi = b.value("phi", 0);
            bin.psi = b.value("psi", 0);
            bin.count = b.value("count", 0);
            bin.rotamers = readRotamers(b.at("rotamers"));
            rr.bins.push_back(std::move(bin));
        }
        lib.residues_.push_back(std::move(rr));
    }
    if (lib.residues_.empty()) throw Error::parse("rotamer pack: `residues` is empty");
    return lib;
}

RotamerLibrary loadRotamerLibrary(const std::filesystem::path& file) {
    std::ifstream in(file);
    if (!in) throw Error::io("cannot open rotamer pack " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse("rotamer pack " + file.string() + ": " + e.what());
    }
    return parseRotamerLibrary(j);
}

const RotamerLibrary* builtinRotamerLibrary() {
    struct Loaded {
        bool ok = false;
        RotamerLibrary lib;
        Loaded() {
            const auto dir = packs::builtinPackDir();
            if (dir.empty()) return;
            const auto file = dir / "rotamers" / "rotamers-pdb-derived-2026.json";
            try {
                lib = loadRotamerLibrary(file);
                ok = true;
            } catch (const Error&) {
            } catch (const std::exception&) {
            }
        }
    };
    static const Loaded loaded;
    return loaded.ok ? &loaded.lib : nullptr;
}

// ---------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------

bool backboneTorsions(const Chain& chain, std::size_t index, double& phi, double& psi) {
    if (index == 0 || index + 1 >= chain.residues.size()) return false;
    const Residue& prev = chain.residues[index - 1];
    const Residue& cur = chain.residues[index];
    const Residue& next = chain.residues[index + 1];
    const Atom* cPrev = prev.atom(" C  ");
    const Atom* n = cur.atom(" N  ");
    const Atom* ca = cur.atom(" CA ");
    const Atom* c = cur.atom(" C  ");
    const Atom* nNext = next.atom(" N  ");
    if (!cPrev || !n || !ca || !c || !nNext) return false;
    const Point3d pC{cPrev->x, cPrev->y, cPrev->z}, pN{n->x, n->y, n->z};
    const Point3d pCa{ca->x, ca->y, ca->z}, pCc{c->x, c->y, c->z};
    const Point3d pNn{nNext->x, nNext->y, nNext->z};
    phi = dihedralDegrees(pC, pN, pCa, pCc);
    psi = dihedralDegrees(pN, pCa, pCc, pNn);
    return true;
}

std::vector<Atom> buildSideChain(const ResidueRotamers& tmpl, const Residue& backbone,
                                 const std::vector<double>& chi) {
    std::vector<Atom> out;
    if (tmpl.build.empty()) return out;

    struct Placed {
        std::string name;
        Point3d p;
    };
    std::vector<Placed> placed;
    for (const char* n : {" N  ", " CA ", " C  ", " O  "}) {
        if (const Atom* a = backbone.atom(n)) placed.push_back({trim(n), {a->x, a->y, a->z}});
    }
    const auto find = [&placed](const std::string& n) -> const Point3d* {
        for (const auto& p : placed)
            if (p.name == n) return &p.p;
        return nullptr;
    };
    if (!find("N") || !find("CA") || !find("C")) return out;   // no backbone, no side chain

    for (const BuildAtom& b : tmpl.build) {
        const Point3d* a0 = find(b.parents[0]);
        const Point3d* a1 = find(b.parents[1]);
        const Point3d* a2 = find(b.parents[2]);
        if (!a0 || !a1 || !a2) break;   // an unbuildable tail is dropped, never guessed
        // chi > 0 means this atom turns with chi_k, and `torsion` is its measured
        // OFFSET from it (0 for the atom that defines chi_k, ~120 for a branch such
        // as LEU CD2). chi == 0 means the torsion is rigid and used as measured.
        double torsion = b.torsion;
        if (b.chi > 0) {
            if (static_cast<std::size_t>(b.chi) > chi.size()) break;
            torsion += chi[static_cast<std::size_t>(b.chi) - 1];
        }
        const Point3d p = placeAtom(*a0, *a1, *a2, b.bondLength, b.bondAngle, torsion);
        placed.push_back({b.atom, p});
        Atom atom;
        // PDB atom names are column-aligned everywhere else in the tree, so the
        // built atoms use the same convention rather than a bare name.
        atom.name = b.atom.size() >= 4 ? b.atom : " " + b.atom + std::string(3 - b.atom.size(), ' ');
        atom.element = b.atom.substr(0, 1);
        atom.x = p.x;
        atom.y = p.y;
        atom.z = p.z;
        out.push_back(std::move(atom));
    }
    return out;
}

int countClashes(const std::vector<Atom>& a, const std::vector<Atom>& b, double overlapFactor) {
    int n = 0;
    for (const Atom& x : a) {
        const double rx = vdwRadius(x.element);
        for (const Atom& y : b) {
            const double cut = overlapFactor * (rx + vdwRadius(y.element));
            const double dx = x.x - y.x, dy = x.y - y.y, dz = x.z - y.z;
            if (dx * dx + dy * dy + dz * dz < cut * cut) ++n;
        }
    }
    return n;
}

// ---------------------------------------------------------------------------
// Dead-end elimination
// ---------------------------------------------------------------------------

double DeeProblem::total(const std::vector<int>& choice) const {
    double t = 0;
    for (std::size_t i = 0; i < choice.size(); ++i) {
        t += self[i][static_cast<std::size_t>(choice[i])];
        for (std::size_t j = i + 1; j < choice.size(); ++j)
            t += pair[i][j][static_cast<std::size_t>(choice[i])][static_cast<std::size_t>(choice[j])];
    }
    return t;
}

DeeResult goldsteinDee(const DeeProblem& problem, std::size_t exhaustiveLimit) {
    DeeResult res;
    const std::size_t np = problem.positions();
    res.alive.resize(np);
    for (std::size_t i = 0; i < np; ++i) res.alive[i].assign(problem.self[i].size(), 1);
    if (np == 0) return res;

    bool changed = true;
    while (changed) {
        changed = false;
        ++res.passes;
        for (std::size_t i = 0; i < np; ++i) {
            const std::size_t nr = problem.self[i].size();
            for (std::size_t r = 0; r < nr; ++r) {
                if (!res.alive[i][r]) continue;
                // Never eliminate the last survivor at a position.
                if (std::count(res.alive[i].begin(), res.alive[i].end(), 1) <= 1) break;
                for (std::size_t t = 0; t < nr; ++t) {
                    if (t == r || !res.alive[i][t]) continue;
                    double bound = problem.self[i][r] - problem.self[i][t];
                    for (std::size_t j = 0; j < np; ++j) {
                        if (j == i) continue;
                        double worst = 0;
                        bool first = true;
                        for (std::size_t s = 0; s < problem.self[j].size(); ++s) {
                            if (!res.alive[j][s]) continue;
                            const double d = problem.pair[i][j][r][s] - problem.pair[i][j][t][s];
                            if (first || d < worst) {
                                worst = d;
                                first = false;
                            }
                        }
                        if (!first) bound += worst;
                    }
                    if (bound > 0) {
                        res.alive[i][r] = 0;
                        ++res.eliminated;
                        changed = true;
                        break;
                    }
                }
            }
        }
    }

    // Final repack over the survivors.
    std::vector<std::vector<int>> survivors(np);
    std::size_t space = 1;
    bool overflow = false;
    for (std::size_t i = 0; i < np; ++i) {
        for (std::size_t r = 0; r < res.alive[i].size(); ++r)
            if (res.alive[i][r]) survivors[i].push_back(static_cast<int>(r));
        if (survivors[i].empty()) survivors[i].push_back(0);
        if (!overflow) {
            if (space > exhaustiveLimit / survivors[i].size()) overflow = true;
            else space *= survivors[i].size();
        }
    }

    res.chosen.assign(np, 0);
    for (std::size_t i = 0; i < np; ++i) res.chosen[i] = survivors[i].front();

    if (!overflow) {
        res.exhaustive = true;
        std::vector<int> idx(np, 0), best;
        double bestScore = 0;
        bool first = true;
        while (true) {
            std::vector<int> choice(np);
            for (std::size_t i = 0; i < np; ++i)
                choice[i] = survivors[i][static_cast<std::size_t>(idx[i])];
            const double s = problem.total(choice);
            if (first || s < bestScore) {
                bestScore = s;
                best = choice;
                first = false;
            }
            std::size_t k = 0;
            for (; k < np; ++k) {
                if (++idx[k] < static_cast<int>(survivors[k].size())) break;
                idx[k] = 0;
            }
            if (k == np) break;
        }
        res.chosen = best;
        res.totalScore = bestScore;
    } else {
        // Iterative single-position descent: with a survivor space this large the
        // exact answer is not affordable, and the result says so via `exhaustive`.
        bool improved = true;
        int guard = 0;
        while (improved && guard++ < 100) {
            improved = false;
            for (std::size_t i = 0; i < np; ++i) {
                int bestR = res.chosen[i];
                double bestS = problem.total(res.chosen);
                for (int r : survivors[i]) {
                    std::vector<int> trial = res.chosen;
                    trial[i] = r;
                    const double s = problem.total(trial);
                    if (s < bestS) {
                        bestS = s;
                        bestR = r;
                        improved = true;
                    }
                }
                res.chosen[i] = bestR;
            }
        }
        res.totalScore = problem.total(res.chosen);
    }
    return res;
}

}  // namespace biocad::bio
