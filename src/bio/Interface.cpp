#include "bio/Interface.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "bio/Imgt.h"
#include "bio/Score.h"

namespace biocad::bio {
namespace {

struct Vec3 {
    double x = 0, y = 0, z = 0;
};

double dist2(const Atom& a, const Atom& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// Atom names are stored VERBATIM with their fixed-column padding (" CA "), because
// the padding encodes the element column. Every name comparison here therefore goes
// through this trim: comparing a.name == "NZ" silently matches nothing.
std::string atomName(const Atom& a) {
    const auto b = a.name.find_first_not_of(' ');
    if (b == std::string::npos) return {};
    const auto e = a.name.find_last_not_of(' ');
    return a.name.substr(b, e - b + 1);
}

bool isHydrogen(const Atom& a) {
    return a.element == "H" || a.element == "D";
}

std::string residueLabel(const std::string& chainId, const Residue& r) {
    std::string s = chainId + ":" + r.name + std::to_string(r.authSeqId);
    if (r.insertionCode != ' ') s.push_back(r.insertionCode);
    return s;
}

std::string residueKey(const std::string& chainId, const Residue& r) {
    std::string s = chainId + ":" + std::to_string(r.authSeqId);
    if (r.insertionCode != ' ') s.push_back(r.insertionCode);
    return s;
}

// A structure containing only the named chains, so SASA of the isolated side is the
// same code path as SASA of the complex - there is one SASA implementation.
Structure subset(const Structure& s, const std::vector<std::string>& chains) {
    Structure out;
    out.id = s.id;
    out.source = s.source;
    const Model* m = s.model(1);
    if (!m) return out;
    Model mm;
    mm.number = m->number;
    for (const auto& c : m->chains)
        if (std::find(chains.begin(), chains.end(), c.id) != chains.end()) mm.chains.push_back(c);
    out.models.push_back(std::move(mm));
    return out;
}

bool aromatic(const std::string& name) {
    return name == "PHE" || name == "TYR" || name == "TRP" || name == "HIS";
}

// Ring centroid of an aromatic side chain, from the atom names the ring is made of.
bool ringCentroid(const Residue& r, Vec3& out) {
    static const std::map<std::string, std::vector<std::string>> rings = {
        {"PHE", {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}},
        {"TYR", {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}},
        {"TRP", {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}},
        {"HIS", {"CG", "ND1", "CD2", "CE1", "NE2"}},
    };
    const auto it = rings.find(r.name);
    if (it == rings.end()) return false;
    int n = 0;
    Vec3 c;
    for (const auto& nm : it->second) {
        const Atom* a = r.atom(nm);
        if (!a) continue;
        c.x += a->x;
        c.y += a->y;
        c.z += a->z;
        ++n;
    }
    if (n < 3) return false;
    out = {c.x / n, c.y / n, c.z / n};
    return true;
}

bool cationCentre(const Residue& r, Vec3& out) {
    if (r.name == "LYS") {
        if (const Atom* a = r.atom("NZ")) {
            out = {a->x, a->y, a->z};
            return true;
        }
        return false;
    }
    if (r.name == "ARG") {
        int n = 0;
        Vec3 c;
        for (const char* nm : {"NE", "CZ", "NH1", "NH2"}) {
            if (const Atom* a = r.atom(nm)) {
                c.x += a->x;
                c.y += a->y;
                c.z += a->z;
                ++n;
            }
        }
        if (n < 2) return false;
        out = {c.x / n, c.y / n, c.z / n};
        return true;
    }
    return false;
}

bool anionAtom(const Residue& r, const Atom& a) {
    const std::string n = atomName(a);
    if (r.name == "ASP") return n == "OD1" || n == "OD2";
    if (r.name == "GLU") return n == "OE1" || n == "OE2";
    return false;
}
bool cationAtom(const Residue& r, const Atom& a) {
    const std::string n = atomName(a);
    if (r.name == "LYS") return n == "NZ";
    if (r.name == "ARG") return n == "NE" || n == "NH1" || n == "NH2";
    if (r.name == "HIS") return n == "ND1" || n == "NE2";
    return false;
}
bool apolar(const std::string& name) {
    static const std::set<std::string> s = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE",
                                            "TRP", "PRO", "CYS", "TYR"};
    return s.count(name) > 0;
}

// Side-chain atoms beyond C-beta, i.e. exactly what an alanine truncation removes.
bool beyondCBeta(const Atom& a) {
    static const std::set<std::string> keep = {"N", "CA", "C", "O", "CB", "OXT"};
    return keep.count(atomName(a)) == 0 && !isHydrogen(a);
}

double totalSasa(const Structure& s) {
    const SasaResult r = sasa(s);
    return r.total.provenance == Provenance::Measured ? r.total.value : 0.0;
}

std::map<std::string, double> relativeSasaMap(const Structure& s) {
    std::map<std::string, double> out;
    const SasaResult r = sasa(s);
    for (const auto& pr : r.perResidue) {
        if (pr.relative.provenance != Provenance::Measured) continue;
        std::string key = pr.key.chainId + ":" + std::to_string(pr.key.authSeqId);
        if (pr.key.insertionCode != ' ') key.push_back(pr.key.insertionCode);
        out[key] = pr.relative.value;
    }
    return out;
}

struct ContactSet {
    std::vector<ResidueContact> contacts;
    std::set<std::string>       residuesA, residuesB;   // labels
    std::set<std::string>       keysA, keysB;
};

ContactSet contactsBetween(const Structure& s, const std::vector<std::string>& a,
                           const std::vector<std::string>& b, const InterfaceOptions& o) {
    ContactSet cs;
    const Model* m = s.model(1);
    if (!m) return cs;
    auto chainsOf = [&](const std::vector<std::string>& ids) {
        std::vector<const Chain*> v;
        for (const auto& c : m->chains)
            if (std::find(ids.begin(), ids.end(), c.id) != ids.end()) v.push_back(&c);
        return v;
    };
    const auto ca = chainsOf(a), cb = chainsOf(b);
    const double cut2 = o.contactCutoff * o.contactCutoff;

    for (const Chain* x : ca) {
        for (const Residue& rx : x->residues) {
            if (rx.oneLetter() == 'X') continue;
            for (const Chain* y : cb) {
                for (const Residue& ry : y->residues) {
                    if (ry.oneLetter() == 'X') continue;
                    ResidueContact rc;
                    double best = 1e18;
                    int pairs = 0;
                    for (const Atom& ax : rx.atoms) {
                        if (isHydrogen(ax)) continue;
                        for (const Atom& ay : ry.atoms) {
                            if (isHydrogen(ay)) continue;
                            const double d2 = dist2(ax, ay);
                            if (d2 > cut2) continue;
                            ++pairs;
                            best = std::min(best, d2);
                            const bool polarX = ax.element == "N" || ax.element == "O";
                            const bool polarY = ay.element == "N" || ay.element == "O";
                            if (polarX && polarY &&
                                d2 <= o.hydrogenBondCutoff * o.hydrogenBondCutoff)
                                rc.hydrogenBond = true;
                            if (d2 <= o.saltBridgeCutoff * o.saltBridgeCutoff &&
                                ((anionAtom(rx, ax) && cationAtom(ry, ay)) ||
                                 (cationAtom(rx, ax) && anionAtom(ry, ay))))
                                rc.saltBridge = true;
                            if (ax.element == "C" && ay.element == "C" && apolar(rx.name) &&
                                apolar(ry.name))
                                rc.hydrophobic = true;
                            if (atomName(ax) == "SG" && atomName(ay) == "SG" &&
                                d2 <= o.disulfideCutoff * o.disulfideCutoff)
                                rc.disulfide = true;
                        }
                    }
                    if (pairs == 0) continue;
                    // pi-pi and cation-pi are centroid geometries and have a longer
                    // reach than the atom cutoff, so they are tested per residue pair.
                    Vec3 gx, gy;
                    if (aromatic(rx.name) && aromatic(ry.name) && ringCentroid(rx, gx) &&
                        ringCentroid(ry, gy)) {
                        const double dx = gx.x - gy.x, dy = gx.y - gy.y, dz = gx.z - gy.z;
                        if (dx * dx + dy * dy + dz * dz <=
                            o.piStackingCutoff * o.piStackingCutoff)
                            rc.piStacking = true;
                    }
                    auto cationPi = [&](const Residue& cat, const Residue& ring) {
                        Vec3 cc, rr;
                        if (!cationCentre(cat, cc) || !ringCentroid(ring, rr)) return false;
                        const double dx = cc.x - rr.x, dy = cc.y - rr.y, dz = cc.z - rr.z;
                        return dx * dx + dy * dy + dz * dz <=
                               o.cationPiCutoff * o.cationPiCutoff;
                    };
                    if (cationPi(rx, ry) || cationPi(ry, rx)) rc.cationPi = true;

                    rc.chainA = x->id;
                    rc.residueA = residueLabel(x->id, rx);
                    rc.chainB = y->id;
                    rc.residueB = residueLabel(y->id, ry);
                    rc.minDistance = std::sqrt(best);
                    rc.atomContacts = pairs;
                    cs.residuesA.insert(rc.residueA);
                    cs.residuesB.insert(rc.residueB);
                    cs.keysA.insert(residueKey(x->id, rx));
                    cs.keysB.insert(residueKey(y->id, ry));
                    cs.contacts.push_back(std::move(rc));
                }
            }
        }
    }
    std::sort(cs.contacts.begin(), cs.contacts.end(),
              [](const ResidueContact& p, const ResidueContact& q) {
                  return p.minDistance < q.minDistance;
              });
    return cs;
}

}  // namespace

std::vector<std::string> parseChainList(const std::string& spec) {
    std::vector<std::string> out;
    if (spec.find(',') != std::string::npos) {
        std::string cur;
        for (char c : spec) {
            if (c == ',') {
                if (!cur.empty()) out.push_back(cur);
                cur.clear();
            } else if (c != ' ') {
                cur.push_back(c);
            }
        }
        if (!cur.empty()) out.push_back(cur);
        return out;
    }
    for (char c : spec)
        if (c != ' ') out.push_back(std::string(1, c));
    return out;
}

InterfaceReport interfaceOf(const Structure& complex, const std::string& chainsA,
                            const std::string& chainsB, const InterfaceOptions& opts) {
    InterfaceReport rep;
    rep.chainsA = chainsA;
    rep.chainsB = chainsB;
    const auto a = parseChainList(chainsA), b = parseChainList(chainsB);
    const Model* m = complex.model(1);
    if (!m) {
        rep.warnings.push_back("the structure has no model 1");
        rep.buriedSurfaceArea = notComputed("a model");
        rep.sasaA = rep.sasaB = rep.sasaComplex = rep.buriedSurfaceArea;
        return rep;
    }
    std::set<std::string> known;
    for (const auto& c : m->chains) known.insert(c.id);
    for (const auto& id : a)
        if (!known.count(id)) rep.warnings.push_back("chain '" + id + "' is not in the model");
    for (const auto& id : b)
        if (!known.count(id)) rep.warnings.push_back("chain '" + id + "' is not in the model");

    const SasaResult probe = sasa(complex);
    rep.sasaParameters = probe.method;

    std::vector<std::string> both = a;
    both.insert(both.end(), b.begin(), b.end());
    const Structure sa = subset(complex, a), sb = subset(complex, b),
                    sab = subset(complex, both);
    const double areaA = totalSasa(sa), areaB = totalSasa(sb), areaAb = totalSasa(sab);
    const std::string src = probe.method;
    rep.sasaA = makeQuantity(areaA, "A^2", 0.0, Provenance::Measured, src);
    rep.sasaB = makeQuantity(areaB, "A^2", 0.0, Provenance::Measured, src);
    rep.sasaComplex = makeQuantity(areaAb, "A^2", 0.0, Provenance::Measured, src);
    rep.buriedSurfaceArea =
        makeQuantity(areaA + areaB - areaAb, "A^2", 0.0, Provenance::Measured,
                     "BSA = SASA(A) + SASA(B) - SASA(AB); " + src +
                         ". This is the TOTAL area buried, i.e. both sides of the interface; "
                         "the per-side value is half of it.");

    const ContactSet cs = contactsBetween(complex, a, b, opts);
    rep.contacts = cs.contacts;

    // Levy's partition needs BOTH the isolated and the complexed accessibility: a
    // residue already buried in its own monomer is support, one that was exposed and
    // becomes buried is core, one that stays exposed is rim.
    const auto relIso = [&] {
        std::map<std::string, double> r = relativeSasaMap(sa);
        for (const auto& [k, v] : relativeSasaMap(sb)) r[k] = v;
        return r;
    }();
    const auto relCx = relativeSasaMap(sab);
    std::set<std::string> interfaceKeys = cs.keysA;
    interfaceKeys.insert(cs.keysB.begin(), cs.keysB.end());
    for (const auto& key : interfaceKeys) {
        const auto iso = relIso.find(key), cx = relCx.find(key);
        if (iso == relIso.end() || cx == relCx.end()) continue;
        if (iso->second < opts.burialThreshold) rep.supportResidues.push_back(key);
        else if (cx->second < opts.burialThreshold) rep.coreResidues.push_back(key);
        else rep.rimResidues.push_back(key);
    }

    // Antibody side: CDR contacts, paratope and epitope. Which side is the antibody
    // is an INPUT; nothing here guesses it from the chain letters.
    if (!opts.antibodyChains.empty()) {
        std::set<std::string> cdrKeys;
        for (const auto& id : opts.antibodyChains) {
            const Chain* ch = nullptr;
            for (const auto& c : m->chains)
                if (c.id == id) ch = &c;
            if (!ch) continue;
            const std::vector<char> sv = sequenceOf(*ch);
            const AbDomain dom = numberDomain(std::string(sv.begin(), sv.end()));
            if (!dom.numbered) {
                rep.warnings.push_back(
                    "chain '" + id +
                    "' was named as an antibody chain but did not number as a V-DOMAIN, so no CDR "
                    "contacts are reported for it");
                continue;
            }
            std::set<int> cdrIdx;
            for (const auto& r : dom.residues)
                if (r.region.rfind("CDR", 0) == 0) cdrIdx.insert(r.sequenceIndex);
            int poly = 0;
            for (const auto& res : ch->residues) {
                if (res.oneLetter() == 'X') continue;
                if (cdrIdx.count(poly)) cdrKeys.insert(residueKey(ch->id, res));
                ++poly;
            }
        }
        const bool aIsAntibody = [&] {
            for (const auto& id : opts.antibodyChains)
                if (std::find(a.begin(), a.end(), id) != a.end()) return true;
            return false;
        }();
        const std::set<std::string>& abKeys = aIsAntibody ? cs.keysA : cs.keysB;
        const std::set<std::string>& agKeys = aIsAntibody ? cs.keysB : cs.keysA;
        for (const auto& k : abKeys) {
            rep.paratope.push_back(k);
            if (cdrKeys.count(k)) rep.cdrContacts.push_back(k);
        }
        for (const auto& k : agKeys) rep.epitope.push_back(k);
        if (rep.paratope.size() > rep.cdrContacts.size())
            rep.warnings.push_back(
                std::to_string(rep.paratope.size() - rep.cdrContacts.size()) +
                " paratope residue(s) are OUTSIDE the CDRs: framework contacts are real and are "
                "reported rather than filtered out.");
    }
    return rep;
}

AlanineScanReport alanineScan(const Structure& complex, const std::string& chainsA,
                              const std::string& chainsB, const InterfaceOptions& opts) {
    AlanineScanReport rep;
    rep.benchmarkName = "none shipped";
    // Never a ddG, and never a correlation this build did not compute.
    rep.benchmarkSpearman = notComputed("a measured benchmark subset");
    rep.disclaimer =
        "This is a GEOMETRIC scan: each side chain is truncated beyond C-beta and the interface is "
        "re-measured. The impact is a unit-free rank ordering of how much interface a side chain "
        "holds, NOT a binding free energy - there is no solvation term, no entropy, no relaxation "
        "of the rest of the structure and no electrostatic model in it. Do not read it as ddG.";
    rep.assumptions.push_back(
        "Truncation keeps N, CA, C, O and CB and deletes every other side-chain atom; Gly and Ala "
        "are skipped because there is nothing beyond CB to remove.");
    rep.assumptions.push_back(
        "The remaining structure is NOT re-minimised, so the lost area is an upper bound on the "
        "area that a real mutant would lose.");

    const InterfaceReport base = interfaceOf(complex, chainsA, chainsB, opts);
    if (base.buriedSurfaceArea.provenance != Provenance::Measured) {
        rep.assumptions.push_back("no interface could be measured, so nothing was scanned");
        return rep;
    }
    const auto a = parseChainList(chainsA), b = parseChainList(chainsB);
    const Model* m = complex.model(1);
    if (!m) return rep;

    // Interface residues only: scanning the whole protein would report thousands of
    // zeroes and bury the answer.
    std::set<std::string> interfaceKeys;
    for (const auto& c : base.contacts) {
        interfaceKeys.insert(c.residueA);
        interfaceKeys.insert(c.residueB);
    }

    int baseHb = 0, baseSb = 0;
    for (const auto& c : base.contacts) {
        if (c.hydrogenBond) ++baseHb;
        if (c.saltBridge) ++baseSb;
    }

    for (const auto& c : m->chains) {
        const bool inA = std::find(a.begin(), a.end(), c.id) != a.end();
        const bool inB = std::find(b.begin(), b.end(), c.id) != b.end();
        if (!inA && !inB) continue;
        for (const auto& res : c.residues) {
            if (res.oneLetter() == 'X') continue;
            if (res.name == "GLY" || res.name == "ALA") continue;
            const std::string label = residueLabel(c.id, res);
            if (!interfaceKeys.count(label)) continue;

            Structure mutant = complex;
            Model* mm = &mutant.models.front();
            for (auto& mc : mm->chains) {
                if (mc.id != c.id) continue;
                for (auto& mr : mc.residues) {
                    if (mr.authSeqId != res.authSeqId || mr.insertionCode != res.insertionCode ||
                        mr.name != res.name)
                        continue;
                    mr.atoms.erase(std::remove_if(mr.atoms.begin(), mr.atoms.end(),
                                                  [](const Atom& at) { return beyondCBeta(at); }),
                                   mr.atoms.end());
                    mr.name = "ALA";
                }
            }
            const InterfaceReport mut = interfaceOf(mutant, chainsA, chainsB, opts);
            int mutHb = 0, mutSb = 0;
            for (const auto& mc : mut.contacts) {
                if (mc.hydrogenBond) ++mutHb;
                if (mc.saltBridge) ++mutSb;
            }
            AlanineScanPosition p;
            p.chain = c.id;
            p.residue = label;
            p.lostBuriedAreaA2 =
                base.buriedSurfaceArea.value - mut.buriedSurfaceArea.value;
            p.lostContacts =
                static_cast<int>(base.contacts.size()) - static_cast<int>(mut.contacts.size());
            p.lostHydrogenBonds = baseHb - mutHb;
            p.lostSaltBridges = baseSb - mutSb;
            rep.positions.push_back(std::move(p));
        }
    }

    // The impact score: a normalised blend of the four losses, deliberately with no
    // unit. makeQuantity() rejects a Heuristic that carries one.
    double maxArea = 0;
    int maxContacts = 0, maxPolar = 0;
    for (const auto& p : rep.positions) {
        maxArea = std::max(maxArea, p.lostBuriedAreaA2);
        maxContacts = std::max(maxContacts, p.lostContacts);
        maxPolar = std::max(maxPolar, p.lostHydrogenBonds + p.lostSaltBridges);
    }
    const std::string src =
        "geometric alanine truncation: 0.5 * lost BSA + 0.3 * lost atom contacts + 0.2 * lost "
        "(H-bonds + salt bridges), each normalised to the largest loss in this scan. Rank ordering "
        "only - unit-free by construction, and NOT a ddG.";
    for (auto& p : rep.positions) {
        // A truncation can UNCOVER area (a neighbouring side chain becomes exposed),
        // which shows up as a negative loss. That is a real geometric outcome and
        // means the side chain holds no interface, so it floors at zero rather than
        // dragging the blended rank below the residues that contribute nothing.
        const double area = maxArea > 0 ? std::max(0.0, p.lostBuriedAreaA2) / maxArea : 0.0;
        const double cont = maxContacts > 0
                                ? std::max(0, p.lostContacts) / static_cast<double>(maxContacts)
                                : 0.0;
        const double pol = maxPolar > 0
                               ? std::max(0, p.lostHydrogenBonds + p.lostSaltBridges) /
                                     static_cast<double>(maxPolar)
                               : 0.0;
        p.impact = makeQuantity(0.5 * area + 0.3 * cont + 0.2 * pol, "", 0.0,
                                Provenance::Heuristic, src);
    }
    std::sort(rep.positions.begin(), rep.positions.end(),
              [](const AlanineScanPosition& p, const AlanineScanPosition& q) {
                  return p.impact.value > q.impact.value;
              });
    return rep;
}

}  // namespace biocad::bio
