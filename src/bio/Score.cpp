#include "bio/Score.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <map>
#include <unordered_map>

namespace biocad::bio {
namespace {

constexpr double kPi = 3.14159265358979323846;

std::string trimName(const std::string& s) {
    const auto first = s.find_first_not_of(' ');
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(' ');
    return s.substr(first, last - first + 1);
}

std::string upper(const std::string& s) {
    std::string out = s;
    for (char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

// Alternate locations are alternative models of the SAME atom. Taking both would double
// every distance and every sphere, so only the blank/primary conformer is used.
bool primaryAltLoc(const Atom& a) { return a.altLoc == ' ' || a.altLoc == 'A'; }

bool isHydrogen(const Atom& a) {
    const std::string e = upper(trimName(a.element));
    return e == "H" || e == "D";
}

const Model* firstModel(const Structure& s) {
    if (s.models.empty()) return nullptr;
    return &s.models.front();
}

struct Vec3 { double x, y, z; };

double dist2(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

}  // namespace

std::string ResidueKey::label() const {
    std::string out = chainId + ":" + std::to_string(authSeqId);
    if (insertionCode != ' ' && insertionCode != '\0') out.push_back(insertionCode);
    return out;
}

// ---------------------------------------------------------------------------
// Pairing
// ---------------------------------------------------------------------------

ResiduePairing pairResidues(const Model& model, const Model& reference) {
    ResiduePairing out;

    // Chain id + author number + insertion code is the only identity both PDB and mmCIF
    // agree on; label numbering does not exist in a PDB file at all.
    struct Key {
        std::string chain;
        int seq;
        char icode;
        bool operator<(const Key& o) const {
            if (chain != o.chain) return chain < o.chain;
            if (seq != o.seq) return seq < o.seq;
            return icode < o.icode;
        }
    };

    std::map<Key, const Residue*> modelIndex;
    for (const Chain& c : model.chains)
        for (const Residue& r : c.residues)
            modelIndex.emplace(Key{c.id, r.authSeqId, r.insertionCode}, &r);

    std::size_t matched = 0;
    for (const Chain& c : reference.chains) {
        for (const Residue& r : c.residues) {
            const Key k{c.id, r.authSeqId, r.insertionCode};
            const auto it = modelIndex.find(k);
            if (it == modelIndex.end()) {
                ++out.unmatchedReference;
                continue;
            }
            out.pairs.push_back(ResiduePair{ResidueKey{c.id, r.authSeqId, r.insertionCode},
                                            it->second, &r});
            ++matched;
        }
    }
    out.unmatchedModel = modelIndex.size() >= matched ? modelIndex.size() - matched : 0;
    return out;
}

// ---------------------------------------------------------------------------
// lDDT
// ---------------------------------------------------------------------------

LddtResult lddt(const Structure& model, const Structure& reference, const LddtOptions& opts) {
    LddtResult out;

    const Model* mm = firstModel(model);
    const Model* rm = firstModel(reference);
    if (mm == nullptr || rm == nullptr || model.atomCount() == 0 || reference.atomCount() == 0) {
        out.global = biocad::notComputed("lDDT needs atoms in both structures; model has "
                                       + std::to_string(model.atomCount()) + ", reference has "
                                       + std::to_string(reference.atomCount()));
        return out;
    }
    if (opts.requireEqualAtomCounts && model.atomCount() != reference.atomCount()) {
        out.global = biocad::notComputed(
            "lDDT compares like with like: model has " + std::to_string(model.atomCount())
            + " atoms, reference has " + std::to_string(reference.atomCount()));
        return out;
    }
    if (opts.tolerances.empty()) {
        out.global = biocad::notComputed("lDDT needs at least one distance tolerance");
        return out;
    }

    const ResiduePairing pairing = pairResidues(*mm, *rm);
    out.unmatchedModelResidues = pairing.unmatchedModel;
    out.unmatchedReferenceResidues = pairing.unmatchedReference;

    // Flat atom lists in matched order. residueOf[i] indexes perResidue.
    std::vector<Vec3> modelXyz, refXyz;
    std::vector<std::size_t> residueOf;
    out.perResidue.reserve(pairing.pairs.size());

    for (const ResiduePair& p : pairing.pairs) {
        const std::size_t ri = out.perResidue.size();
        out.perResidue.push_back(ResidueLddt{p.key, 0, 0.0});

        std::unordered_map<std::string, const Atom*> byName;
        for (const Atom& a : p.model->atoms) {
            if (!primaryAltLoc(a)) continue;
            byName.emplace(trimName(a.name), &a);
        }
        for (const Atom& ra : p.reference->atoms) {
            if (!primaryAltLoc(ra)) continue;
            if (!opts.includeHydrogens && isHydrogen(ra)) continue;
            if (!opts.includeHetatm && ra.hetatm) continue;
            const auto it = byName.find(trimName(ra.name));
            if (it == byName.end()) continue;
            const Atom& ma = *it->second;
            refXyz.push_back(Vec3{ra.x, ra.y, ra.z});
            modelXyz.push_back(Vec3{ma.x, ma.y, ma.z});
            residueOf.push_back(ri);
        }
    }

    if (modelXyz.size() < 2) {
        out.global = biocad::notComputed(
            "lDDT matched only " + std::to_string(modelXyz.size())
            + " atoms by (chain, authSeqId, insertionCode, atom name)");
        return out;
    }

    const double r0sq = opts.inclusionRadius * opts.inclusionRadius;
    const std::size_t nTol = opts.tolerances.size();

    std::vector<double> residuePreserved(out.perResidue.size(), 0.0);
    std::vector<double> residueTotal(out.perResidue.size(), 0.0);
    double preserved = 0.0;
    double total = 0.0;

    const std::size_t n = modelXyz.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            // Intra-residue distances are fixed by chemistry, not by modelling quality.
            if (residueOf[i] == residueOf[j]) continue;
            const double dr2 = dist2(refXyz[i], refXyz[j]);
            if (dr2 > r0sq) continue;
            const double dr = std::sqrt(dr2);
            const double dm = std::sqrt(dist2(modelXyz[i], modelXyz[j]));
            const double delta = std::fabs(dm - dr);

            std::size_t hits = 0;
            for (double tol : opts.tolerances)
                if (delta < tol) ++hits;

            ++out.consideredPairs;
            const double h = static_cast<double>(hits);
            const double t = static_cast<double>(nTol);
            preserved += h;
            total += t;
            residuePreserved[residueOf[i]] += h;
            residueTotal[residueOf[i]] += t;
            residuePreserved[residueOf[j]] += h;
            residueTotal[residueOf[j]] += t;
            ++out.perResidue[residueOf[i]].distancePairs;
            ++out.perResidue[residueOf[j]].distancePairs;
        }
    }

    if (total == 0.0) {
        out.global = biocad::notComputed("lDDT found no reference atom pairs within "
                                       + std::to_string(opts.inclusionRadius) + " A");
        return out;
    }

    for (std::size_t i = 0; i < out.perResidue.size(); ++i)
        out.perResidue[i].score =
            residueTotal[i] > 0.0 ? residuePreserved[i] / residueTotal[i] : 0.0;

    std::string src = "lDDT (Mariani 2013), superposition-free, R0="
                      + std::to_string(opts.inclusionRadius) + " A, tolerances";
    for (double tol : opts.tolerances) src += " " + std::to_string(tol);
    src += " A, " + std::to_string(out.consideredPairs) + " reference pairs, hydrogens "
           + (opts.includeHydrogens ? "included" : "ignored");

    // Unitless fraction, and exact over the given coordinates - never a prediction.
    out.global = biocad::makeQuantity(preserved / total, "", 0.0, biocad::Provenance::Measured, src);
    return out;
}

// ---------------------------------------------------------------------------
// SASA
// ---------------------------------------------------------------------------

double vdwRadius(const std::string& element) {
    const std::string e = upper(trimName(element));
    // Bondi 1964 / Rowland & Taylor 1996 values, the set every structural tool quotes.
    if (e == "C") return 1.70;
    if (e == "N") return 1.55;
    if (e == "O") return 1.52;
    if (e == "S") return 1.80;
    if (e == "P") return 1.80;
    if (e == "SE") return 1.90;
    if (e == "H" || e == "D") return 1.20;
    // Documented fallback: an unknown element is given the sulfur radius rather than zero,
    // so it still occludes its neighbours instead of silently opening a hole in the surface.
    return 1.80;
}

double maxAccessibility(const std::string& residueName) {
    // Tien et al. 2013, PLoS ONE 8:e80635, "theoretical" Gly-X-Gly maxima, A^2.
    static const std::map<std::string, double> kMax{
        {"ALA", 129.0}, {"ARG", 274.0}, {"ASN", 195.0}, {"ASP", 193.0}, {"CYS", 167.0},
        {"GLU", 223.0}, {"GLN", 225.0}, {"GLY", 104.0}, {"HIS", 224.0}, {"ILE", 197.0},
        {"LEU", 201.0}, {"LYS", 236.0}, {"MET", 224.0}, {"PHE", 240.0}, {"PRO", 159.0},
        {"SER", 155.0}, {"THR", 172.0}, {"TRP", 285.0}, {"TYR", 263.0}, {"VAL", 174.0},
    };
    const auto it = kMax.find(upper(trimName(residueName)));
    return it == kMax.end() ? 0.0 : it->second;
}

namespace {

// Golden-spiral (Fibonacci) points on the unit sphere: near-uniform for any N, unlike a
// lat/long grid which clusters at the poles and biases small spheres.
std::vector<Vec3> spherePoints(int n) {
    std::vector<Vec3> pts;
    pts.reserve(static_cast<std::size_t>(n));
    const double golden = kPi * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < n; ++i) {
        const double y = 1.0 - 2.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double phi = static_cast<double>(i) * golden;
        pts.push_back(Vec3{std::cos(phi) * r, y, std::sin(phi) * r});
    }
    return pts;
}

std::string formatDouble(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    return buf;
}

}  // namespace

SasaResult sasa(const Structure& s, const SasaOptions& opts) {
    SasaResult out;
    // Every SASA number in this program carries this string. Without it a SASA value is not
    // reproducible: probe radius, point count, radii set and hydrogen policy each move it.
    out.method = "Shrake-Rupley, probe " + formatDouble(opts.probeRadius, 2) + " A, "
                 + std::to_string(opts.points) + " golden-spiral points, radii Bondi "
                 + "(C 1.70, N 1.55, O 1.52, S 1.80, P 1.80, SE 1.90, H 1.20, other 1.80), "
                 + "hydrogens " + (opts.includeHydrogens ? "included (all-atom)"
                                                         : "ignored (united-atom)");

    const Model* m = firstModel(s);
    if (m == nullptr || s.atomCount() == 0) {
        out.total = biocad::notComputed("SASA needs at least one atom; the structure has none");
        return out;
    }
    if (opts.points < 1) {
        out.total = biocad::notComputed("SASA needs at least one quadrature point");
        return out;
    }

    struct Sphere { Vec3 c; double r; std::size_t residue; };
    std::vector<Sphere> spheres;

    for (const Chain& c : m->chains) {
        for (const Residue& r : c.residues) {
            const std::size_t ri = out.perResidue.size();
            out.perResidue.push_back(ResidueSasa{
                ResidueKey{c.id, r.authSeqId, r.insertionCode}, r.name,
                biocad::makeQuantity(0.0, "A^2", 0.0, biocad::Provenance::Measured, out.method),
                biocad::notComputed("no Tien 2013 maximum accessibility for residue " + r.name)});
            for (const Atom& a : r.atoms) {
                if (!primaryAltLoc(a)) continue;
                if (!opts.includeHydrogens && isHydrogen(a)) continue;
                if (!opts.includeHetatm && a.hetatm) continue;
                spheres.push_back(Sphere{Vec3{a.x, a.y, a.z},
                                         vdwRadius(a.element) + opts.probeRadius, ri});
            }
        }
    }

    if (spheres.empty()) {
        out.total = biocad::notComputed(
            "SASA found no eligible atoms (hydrogens "
            + std::string(opts.includeHydrogens ? "included" : "ignored") + ", hetatm "
            + std::string(opts.includeHetatm ? "included" : "ignored") + ")");
        return out;
    }

    const std::vector<Vec3> unit = spherePoints(opts.points);
    const double invP = 1.0 / static_cast<double>(opts.points);

    // Neighbour prefilter: only spheres whose expanded radii overlap can occlude each other.
    // Brute force over the O(n^2) pair test is one squared distance per pair and is far
    // cheaper than the O(n * points * neighbours) point test it prunes.
    const std::size_t n = spheres.size();
    std::vector<std::vector<std::size_t>> neighbours(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double cut = spheres[i].r + spheres[j].r;
            if (dist2(spheres[i].c, spheres[j].c) < cut * cut) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }
        }
    }

    double totalArea = 0.0;
    std::vector<double> residueArea(out.perResidue.size(), 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const Sphere& si = spheres[i];
        int accessible = 0;
        for (const Vec3& u : unit) {
            const Vec3 p{si.c.x + si.r * u.x, si.c.y + si.r * u.y, si.c.z + si.r * u.z};
            bool buried = false;
            for (std::size_t j : neighbours[i]) {
                const Sphere& sj = spheres[j];
                if (dist2(p, sj.c) < sj.r * sj.r) { buried = true; break; }
            }
            if (!buried) ++accessible;
        }
        // SASA_i = 4 pi R_i^2 * A_i / P, with R_i the vdW radius PLUS the probe radius:
        // the accessible surface is the locus of the probe centre, not of the atom skin.
        const double area = 4.0 * kPi * si.r * si.r * static_cast<double>(accessible) * invP;
        totalArea += area;
        residueArea[si.residue] += area;
    }

    for (std::size_t i = 0; i < out.perResidue.size(); ++i) {
        ResidueSasa& rs = out.perResidue[i];
        rs.absolute = biocad::makeQuantity(residueArea[i], "A^2", 0.0, biocad::Provenance::Measured,
                                         out.method);
        const double maxAcc = maxAccessibility(rs.residueName);
        if (maxAcc > 0.0) {
            rs.relative = biocad::makeQuantity(residueArea[i] / maxAcc, "", 0.0,
                                             biocad::Provenance::Measured,
                                             out.method + "; relative to Tien 2013 theoretical "
                                                          "maximum " + formatDouble(maxAcc, 1)
                                                 + " A^2");
        }
    }

    out.total = biocad::makeQuantity(totalArea, "A^2", 0.0, biocad::Provenance::Measured, out.method);
    return out;
}

}  // namespace biocad::bio
