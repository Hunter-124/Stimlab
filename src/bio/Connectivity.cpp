#include "bio/Connectivity.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>

#include <nlohmann/json.hpp>

#include "core/Assets.h"

namespace biocad::bio {
namespace {

using nlohmann::json;

std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string trimUpper(const std::string& s) {
    const auto b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t");
    std::string out = s.substr(b, e - b + 1);
    for (char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return out;
}

BondOrder parseOrder(const std::string& s) {
    if (s == "2") return BondOrder::Double;
    if (s == "ar") return BondOrder::Aromatic;
    return BondOrder::Single;
}

// Two atoms may be bonded only when their alternate-location indicators are compatible: a
// blank altLoc belongs to every conformer, and two lettered altLocs must be the SAME letter.
// Without this rule the A and B copies of a disordered side chain get cross-bonded and the
// residue renders with an extra branch that exists in no conformer.
bool altLocCompatible(char a, char b) { return a == ' ' || b == ' ' || a == b; }

double distance(const Atom& a, const Atom& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

ResidueTemplatePack loadPack(const std::filesystem::path& path) {
    ResidueTemplatePack pack;
    if (path.empty()) {
        pack.errors.push_back("no asset tree found: assets/packs/structure is unavailable");
        return pack;
    }
    const std::string text = readFile(path);
    if (text.empty()) {
        pack.errors.push_back("could not read " + path.string());
        return pack;
    }
    json j;
    try {
        j = json::parse(text);
    } catch (const json::exception& e) {
        pack.errors.push_back(std::string("residue-templates.json is not valid JSON: ") +
                              e.what());
        return pack;
    }
    pack.schemaVersion = j.value("schemaVersion", 0);
    if (pack.schemaVersion != 1) {
        pack.errors.push_back("residue-templates.json schemaVersion " +
                              std::to_string(pack.schemaVersion) + " is not supported (expected 1)");
        return pack;
    }
    try {
        for (const auto& r : j.at("residues")) {
            ResidueTemplate t;
            t.name = trimUpper(r.at("name").get<std::string>());
            t.kind = r.value("kind", "protein");
            for (const auto& b : r.at("bonds")) {
                TemplateBond tb;
                tb.a = trimUpper(b.at(0).get<std::string>());
                tb.b = trimUpper(b.at(1).get<std::string>());
                tb.order = parseOrder(b.at(2).get<std::string>());
                t.bonds.push_back(std::move(tb));
            }
            pack.residues.push_back(std::move(t));
        }
        for (const auto& [k, v] : j.at("aliases").items())
            pack.aliases.emplace_back(trimUpper(k), trimUpper(v.get<std::string>()));
        auto rule = [&](const char* key) {
            const auto& n = j.at("linkage").at(key);
            LinkageRule lr;
            lr.from = trimUpper(n.at("from").get<std::string>());
            lr.to = trimUpper(n.at("to").get<std::string>());
            lr.minAngstrom = n.at("minAngstrom").get<double>();
            lr.maxAngstrom = n.at("maxAngstrom").get<double>();
            return lr;
        };
        pack.protein = rule("protein");
        pack.nucleic = rule("nucleic");
        pack.disulfide = rule("disulfide");
        for (const auto& s : j.at("solvent")) pack.solvent.push_back(trimUpper(s.get<std::string>()));
    } catch (const json::exception& e) {
        pack.errors.push_back(std::string("residue-templates.json is missing a field: ") + e.what());
        return pack;
    }
    pack.ok = !pack.residues.empty();
    return pack;
}

// Atom names indexed by trimmed uppercase name. A name may repeat inside one residue when the
// side chain is disordered (OG with altLoc A and B), so the value is a list.
using NameIndex = std::unordered_map<std::string, std::vector<int>>;

NameIndex indexAtoms(const Residue& r) {
    NameIndex idx;
    for (int i = 0; i < static_cast<int>(r.atoms.size()); ++i)
        idx[trimUpper(r.atoms[i].name)].push_back(i);
    return idx;
}

// One representative atom for an inter-residue link: the blank-altLoc copy when there is one,
// otherwise the first. Linking every altLoc pair would multiply the backbone.
const int* representative(const NameIndex& idx, const Residue& r, const std::string& name) {
    const auto it = idx.find(name);
    if (it == idx.end() || it->second.empty()) return nullptr;
    for (const int& i : it->second)
        if (r.atoms[static_cast<std::size_t>(i)].altLoc == ' ') return &i;
    return &it->second.front();
}

}  // namespace

const ResidueTemplate* ResidueTemplatePack::find(const std::string& residueName) const {
    const std::string key = trimUpper(residueName);
    std::string target = key;
    for (const auto& [from, to] : aliases) {
        if (from == key) {
            target = to;
            break;
        }
    }
    for (const ResidueTemplate& t : residues)
        if (t.name == target) return &t;
    return nullptr;
}

bool ResidueTemplatePack::isSolvent(const std::string& residueName) const {
    const std::string key = trimUpper(residueName);
    return std::find(solvent.begin(), solvent.end(), key) != solvent.end();
}

ResidueTemplatePack loadResidueTemplates(const std::filesystem::path& file) {
    return loadPack(file);
}

const ResidueTemplatePack& residueTemplates() {
    const std::filesystem::path dir = core::assetDir("packs/structure");
    static const ResidueTemplatePack pack =
        loadPack(dir.empty() ? std::filesystem::path{} : dir / "residue-templates.json");
    return pack;
}

ConnectivityResult connect(const Model& model, const ConnectivityOptions& options) {
    ConnectivityResult out;
    const ResidueTemplatePack& pack = residueTemplates();
    for (const std::string& e : pack.errors) out.diagnostics.warnings.push_back(e);

    // degree[chain][residue][atom]; parallel to the hierarchy so a degree-0 atom can be
    // reported with the residue it belongs to.
    std::vector<std::vector<std::vector<int>>> degree(model.chains.size());
    std::vector<std::vector<NameIndex>> names(model.chains.size());
    std::unordered_map<std::string, UnknownResidue> unknown;

    for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
        const Chain& chain = model.chains[ci];
        degree[ci].resize(chain.residues.size());
        names[ci].reserve(chain.residues.size());
        for (std::size_t ri = 0; ri < chain.residues.size(); ++ri) {
            const Residue& res = chain.residues[ri];
            degree[ci][ri].assign(res.atoms.size(), 0);
            names[ci].push_back(indexAtoms(res));
        }
    }

    auto addBond = [&](std::size_t ci, std::size_t ri, int ai, std::size_t cj, std::size_t rj,
                       int aj, BondOrder order, BondKind kind) {
        StructureBond b;
        b.a = {static_cast<int>(ci), static_cast<int>(ri), ai};
        b.b = {static_cast<int>(cj), static_cast<int>(rj), aj};
        b.order = order;
        b.kind = kind;
        out.bonds.push_back(b);
        ++degree[ci][ri][static_cast<std::size_t>(ai)];
        ++degree[cj][rj][static_cast<std::size_t>(aj)];
    };

    // --- Intra-residue: templates by atom name ---
    for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
        const Chain& chain = model.chains[ci];
        for (std::size_t ri = 0; ri < chain.residues.size(); ++ri) {
            const Residue& res = chain.residues[ri];
            const ResidueTemplate* t = pack.find(res.name);
            if (!t) {
                if (pack.isSolvent(res.name)) {
                    ++out.diagnostics.solventResidues;
                } else if (res.atoms.size() <= 1) {
                    // A single-atom residue is an ion or a lone solvent oxygen: there is no bond
                    // to make, so calling it "unknown" would be a warning about nothing.
                    ++out.diagnostics.monatomicResidues;
                } else {
                    UnknownResidue& u = unknown[trimUpper(res.name)];
                    u.name = trimUpper(res.name);
                    ++u.count;
                    u.atoms += static_cast<int>(res.atoms.size());
                }
                continue;
            }
            const NameIndex& idx = names[ci][ri];
            for (const TemplateBond& tb : t->bonds) {
                const auto ita = idx.find(tb.a);
                const auto itb = idx.find(tb.b);
                if (ita == idx.end() || itb == idx.end()) {
                    // OXT/OP3 are terminal-only, and a disordered side chain simply stops. Both
                    // are "the atom is not in the file", which is counted, never invented.
                    ++out.diagnostics.missingTemplateAtoms;
                    continue;
                }
                for (int ai : ita->second) {
                    for (int bi : itb->second) {
                        if (!altLocCompatible(res.atoms[static_cast<std::size_t>(ai)].altLoc,
                                              res.atoms[static_cast<std::size_t>(bi)].altLoc))
                            continue;
                        addBond(ci, ri, ai, ci, ri, bi, tb.order, BondKind::Template);
                        ++out.diagnostics.templateBonds;
                    }
                }
            }
        }
    }

    // --- Inter-residue polymer links, and the gaps where there is no link ---
    for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
        const Chain& chain = model.chains[ci];
        for (std::size_t ri = 0; ri + 1 < chain.residues.size(); ++ri) {
            const Residue& a = chain.residues[ri];
            const Residue& b = chain.residues[ri + 1];
            const ResidueTemplate* ta = pack.find(a.name);
            const ResidueTemplate* tb = pack.find(b.name);
            if (!ta || !tb || ta->kind != tb->kind) continue;   // not a polymer step
            const bool protein = ta->kind == "protein";
            const LinkageRule& rule = protein ? pack.protein : pack.nucleic;

            const int* ia = representative(names[ci][ri], a, rule.from);
            const int* ib = representative(names[ci][ri + 1], b, rule.to);
            ChainGap gap;
            gap.chainId = chain.id;
            gap.fromSeqId = a.authSeqId;
            gap.fromInsertionCode = a.insertionCode;
            gap.toSeqId = b.authSeqId;
            gap.toInsertionCode = b.insertionCode;
            if (!ia || !ib) {
                gap.reason = "linking atom " + (ia ? rule.to : rule.from) + " is absent";
                out.diagnostics.gaps.push_back(gap);
                continue;
            }
            const double d = distance(a.atoms[static_cast<std::size_t>(*ia)],
                                      b.atoms[static_cast<std::size_t>(*ib)]);
            if (d < rule.minAngstrom || d > rule.maxAngstrom) {
                gap.distance = d;
                gap.reason = rule.from + ".." + rule.to + " is " +
                             std::to_string(d) + " A, outside the bonding window";
                out.diagnostics.gaps.push_back(gap);
                continue;
            }
            addBond(ci, ri, *ia, ci, ri + 1, *ib, BondOrder::Single,
                    protein ? BondKind::PeptideLink : BondKind::NucleicLink);
            ++out.diagnostics.linkBonds;
        }
    }

    // --- Disulfides, by distance, across chains ---
    // SSBOND / _struct_conn is optional and routinely absent, so the geometry is the source.
    // It is reported as its own BondKind so it is never confused with a template bond.
    if (options.findDisulfides) {
        struct Sg { std::size_t ci, ri; int ai; };
        std::vector<Sg> sgs;
        for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
            const Chain& chain = model.chains[ci];
            for (std::size_t ri = 0; ri < chain.residues.size(); ++ri) {
                const ResidueTemplate* t = pack.find(chain.residues[ri].name);
                if (!t || t->name != "CYS") continue;
                const auto it = names[ci][ri].find(pack.disulfide.from);
                if (it == names[ci][ri].end()) continue;
                for (int ai : it->second) sgs.push_back({ci, ri, ai});
            }
        }
        for (std::size_t i = 0; i < sgs.size(); ++i) {
            for (std::size_t j = i + 1; j < sgs.size(); ++j) {
                if (sgs[i].ci == sgs[j].ci && sgs[i].ri == sgs[j].ri) continue;
                const Atom& x = model.chains[sgs[i].ci]
                                    .residues[sgs[i].ri]
                                    .atoms[static_cast<std::size_t>(sgs[i].ai)];
                const Atom& y = model.chains[sgs[j].ci]
                                    .residues[sgs[j].ri]
                                    .atoms[static_cast<std::size_t>(sgs[j].ai)];
                if (!altLocCompatible(x.altLoc, y.altLoc)) continue;
                const double d = distance(x, y);
                if (d < pack.disulfide.minAngstrom || d > pack.disulfide.maxAngstrom) continue;
                addBond(sgs[i].ci, sgs[i].ri, sgs[i].ai, sgs[j].ci, sgs[j].ri, sgs[j].ai,
                        BondOrder::Single, BondKind::Disulfide);
                ++out.diagnostics.disulfides;
            }
        }
    }

    // --- Opt-in distance fallback, for residues with no template only ---
    if (options.inferByDistance) {
        const double cut = options.inferMaxAngstrom;
        for (std::size_t ci = 0; ci < model.chains.size(); ++ci) {
            const Chain& chain = model.chains[ci];
            for (std::size_t ri = 0; ri < chain.residues.size(); ++ri) {
                const Residue& res = chain.residues[ri];
                if (pack.find(res.name) || pack.isSolvent(res.name) || res.atoms.size() <= 1)
                    continue;
                for (std::size_t i = 0; i < res.atoms.size(); ++i) {
                    for (std::size_t j = i + 1; j < res.atoms.size(); ++j) {
                        if (!altLocCompatible(res.atoms[i].altLoc, res.atoms[j].altLoc)) continue;
                        const double d = distance(res.atoms[i], res.atoms[j]);
                        if (d < 0.4 || d > cut) continue;
                        addBond(ci, ri, static_cast<int>(i), ci, ri, static_cast<int>(j),
                                BondOrder::Single, BondKind::DistanceInferred);
                        ++out.diagnostics.inferredBonds;
                    }
                }
            }
        }
    }

    // --- Diagnostics roll-up ---
    for (auto& [key, u] : unknown) {
        (void)key;
        out.diagnostics.unknownResidues.push_back(u);
    }
    std::sort(out.diagnostics.unknownResidues.begin(), out.diagnostics.unknownResidues.end(),
              [](const UnknownResidue& a, const UnknownResidue& b) { return a.name < b.name; });

    for (std::size_t ci = 0; ci < model.chains.size(); ++ci)
        for (std::size_t ri = 0; ri < model.chains[ci].residues.size(); ++ri)
            for (int d : degree[ci][ri])
                if (d == 0) ++out.diagnostics.unbondedAtoms;

    for (const UnknownResidue& u : out.diagnostics.unknownResidues) {
        out.diagnostics.warnings.push_back(
            "no connectivity template for residue " + u.name + " (" + std::to_string(u.count) +
            " occurrence(s), " + std::to_string(u.atoms) +
            " atoms): left unbonded. Enable the distance fallback to bond it as a guess.");
    }
    for (const ChainGap& g : out.diagnostics.gaps) {
        out.diagnostics.warnings.push_back("chain " + g.chainId + " break between " +
                                           std::to_string(g.fromSeqId) + g.fromInsertionCode +
                                           " and " + std::to_string(g.toSeqId) +
                                           g.toInsertionCode + ": " + g.reason);
    }
    return out;
}

}  // namespace biocad::bio
