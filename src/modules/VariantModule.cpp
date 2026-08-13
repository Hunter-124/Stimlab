#include "modules/VariantModule.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "bio/Conservation.h"
#include "bio/Rotamer.h"
#include "bio/Score.h"
#include "modules/BioModules.h"

namespace biocad {
namespace {

const char* kSiftSource =
    "SIFT-style tolerance index p(mutant)/max_y p(y) over the supplied homolog set "
    "(Ng & Henikoff 2003); deleterious below 0.05";
const char* kProveanSource =
    "PROVEAN-style delta alignment score, BLOSUM62 (Choi et al. 2012, PLoS ONE "
    "7:e46688; threshold -2.282 at 79.05% balanced accuracy on 58,684 human variants)";

std::string shortfall(const ConservationProfile& p) {
    return "at least " + std::to_string(p.minimumHomologsRequired) + " homologs (" +
           std::to_string(p.homologs.sequenceCount) + " supplied)";
}

// Heavy atoms of a residue's backbone, which no rotamer choice can move and which
// are therefore excluded from the clash count of the side chains attached to or
// adjacent to them.
bool isBackboneName(const std::string& n) {
    const std::string t = n.substr(n.find_first_not_of(' ') == std::string::npos
                                       ? 0
                                       : n.find_first_not_of(' '));
    const std::string s = t.substr(0, t.find_last_not_of(' ') + 1);
    return s == "N" || s == "CA" || s == "C" || s == "O" || s == "OXT";
}

const bio::Atom* anchorAtom(const bio::Residue& r) {
    if (const bio::Atom* cb = r.atom(" CB ")) return cb;
    return r.atom(" CA ");
}

struct Site {
    const bio::Chain* chain = nullptr;
    std::size_t index = 0;
    std::string label;              // "A:112"
    std::string residueName;        // the type being packed at this site
    std::vector<bio::RotamerEntry> rotamers;
    bool backboneIndependent = true;
    std::vector<std::vector<bio::Atom>> built;   // one atom set per rotamer
};

}  // namespace

std::string threeLetterOf(char oneLetter) {
    static const char* kOne = "ACDEFGHIKLMNPQRSTVWY";
    static const char* kThree[] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS",
                                   "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
                                   "ARG", "SER", "THR", "VAL", "TRP", "TYR"};
    if (oneLetter >= 'a' && oneLetter <= 'z')
        oneLetter = static_cast<char>(oneLetter - 'a' + 'A');
    for (int i = 0; kOne[i]; ++i)
        if (kOne[i] == oneLetter) return kThree[i];
    return "UNK";
}

ConservationProfile RealVariants::conservation(const std::string& query,
                                               const std::vector<std::string>& homologs) const {
    const bio::SubstitutionMatrix* m = bioadapt::blosum62();
    if (!m) {
        ConservationProfile p;
        p.query = query;
        p.minimumHomologsRequired = bio::kMinimumHomologs;
        p.homologs.sequenceCount = static_cast<int>(homologs.size());
        p.usable = false;
        p.warnings.push_back(
            "assets/packs/matrices/blosum62.json could not be loaded, so no alignment was "
            "attempted and no conservation was computed. Nothing is substituted for it.");
        return p;
    }
    bio::ProfileOptions opts;
    opts.source = "user-supplied homolog set (" + std::to_string(homologs.size()) +
                  " sequences), aligned with " + bioadapt::kAlignSource();
    return bio::buildProfile("query", query, homologs, *m, opts);
}

VariantScore RealVariants::score(const ConservationProfile& profile, int position,
                                 char mutant) const {
    VariantScore out;
    out.position = position;
    out.mutant = mutant;
    out.homologs = profile.homologs;

    const bio::SubstitutionMatrix* m = bioadapt::blosum62();
    const ConservationColumn* col = nullptr;
    for (const auto& c : profile.columns)
        if (c.position == position) col = &c;
    if (col) {
        out.wildType = col->queryResidue;
    } else if (position >= 1 && static_cast<std::size_t>(position) <= profile.query.size()) {
        // An unusable profile has no columns, but the query sequence is still known,
        // and the BLOSUM62 delta does not depend on the homolog set. Refusing it
        // here would throw away the cheap honest number along with the expensive
        // uncertain one.
        out.wildType = profile.query[static_cast<std::size_t>(position) - 1];
    }

    if (!m) {
        const std::string missing = "assets/packs/matrices/blosum62.json";
        out.blosum62Delta = notComputed(missing);
        out.siftScore = notComputed(missing);
        out.proveanScore = notComputed(missing);
        out.columnEntropy = notComputed(missing);
        out.interpretation = "The BLOSUM62 matrix could not be loaded, so nothing was scored.";
        return out;
    }

    // The BLOSUM62 delta does not depend on the homolog set: it is score(wt, mut)
    // minus score(wt, wt) straight out of the matrix, so it is reported even when
    // the conservation half is refused.
    if (out.wildType != 'X') {
        out.blosum62Delta = makeQuantity(
            static_cast<double>(bio::blosum62Delta(*m, out.wildType, mutant)), "half-bits", 0.0,
            Provenance::Measured,
            "BLOSUM62 table lookup: score(" + std::string(1, out.wildType) + "," +
                std::string(1, mutant) + ") - score(" + std::string(1, out.wildType) + "," +
                std::string(1, out.wildType) + "), " + m->source());
    } else {
        out.blosum62Delta = notComputed("a wild-type residue at position " +
                                        std::to_string(position));
    }

    if (!profile.usable || !col) {
        // Order matters: an unusable profile has no columns at all, so testing the
        // missing column first would report "position N is not in a 0-column
        // profile" and hide the reason that actually applies.
        const std::string missing =
            !profile.usable
                ? shortfall(profile)
                : ("position " + std::to_string(position) + " inside the profile (it has " +
                   std::to_string(profile.columns.size()) + " columns)");
        out.siftScore = notComputed(missing);
        out.proveanScore = notComputed(missing);
        out.columnEntropy = notComputed(missing);
        out.warnings = profile.warnings;
        out.interpretation =
            "No conservation score was produced because the profile is not usable: " + missing +
            ". The BLOSUM62 delta above is a matrix lookup and stands on its own; it is not a "
            "conservation measurement and it is not a pathogenicity call.";
        return out;
    }

    const double sift = bio::siftScore(*col, mutant);
    const double provean = bio::proveanDelta(*col, out.wildType, mutant, *m);
    // Both are dimensionless ratios/score differences over a user-supplied set:
    // no published model ran, so neither is Predicted; they rank, so Heuristic
    // with an empty unit, which makeQuantity enforces.
    out.siftScore = makeQuantity(sift, "", 0.0, Provenance::Heuristic,
                                 std::string(kSiftSource) + ", " +
                                     std::to_string(profile.homologs.sequenceCount) +
                                     " homologs, median identity " +
                                     std::to_string(static_cast<int>(std::lround(
                                         profile.homologs.medianIdentityPct))) +
                                     "%");
    out.proveanScore = makeQuantity(provean, "", 0.0, Provenance::Heuristic, kProveanSource);
    out.columnEntropy = makeQuantity(col->shannonEntropy, "bits", 0.0, Provenance::Measured,
                                     "Shannon entropy of the observed weighted column "
                                     "distribution over the 20 standard residues; 0 is fully "
                                     "conserved, log2(20) = 4.3219 is uniform");
    out.warnings = profile.warnings;
    out.interpretation =
        std::string("SIFT-style tolerance index ") + std::to_string(sift) +
        " (deleterious below 0.05) and PROVEAN-style delta " + std::to_string(provean) +
        " (deleterious below -2.282), both computed from " +
        std::to_string(profile.homologs.sequenceCount) + " user-supplied homologs at median " +
        std::to_string(static_cast<int>(std::lround(profile.homologs.medianIdentityPct))) +
        "% identity. These describe how unusual this residue would be in THIS alignment. "
        "They are not a pathogenicity call, they are not a clinical interpretation, and a "
        "score either side of its threshold is not a verdict about a person.";
    return out;
}

RotamerRebuild RealVariants::rebuild(const bio::Structure& s, const std::string& chainId,
                                     int residueNumber, char mutant) const {
    RotamerRebuild out;
    out.position = residueNumber;
    out.mutant = mutant;
    out.provenance = Provenance::Model;

    const bio::RotamerLibrary* lib = bio::builtinRotamerLibrary();
    if (!lib) {
        out.warnings.push_back(
            "assets/packs/rotamers/rotamers-pdb-derived-2026.json could not be loaded, so no "
            "side chain was built. No rotamer angle is guessed in its place.");
        return out;
    }
    out.rotamerLibrarySource = lib->name() + " (" + std::to_string(lib->datasetEntryCount()) +
                               " PDB entries, " + std::to_string(lib->residuesMeasured()) +
                               " side chains measured)";

    const bio::Model* model = s.model(1);
    if (!model) {
        out.warnings.push_back("The structure has no model 1, so there was nothing to rebuild.");
        return out;
    }
    const bio::Chain* target = nullptr;
    std::size_t targetIndex = 0;
    for (const auto& c : model->chains) {
        if (c.id != chainId) continue;
        for (std::size_t i = 0; i < c.residues.size(); ++i) {
            if (c.residues[i].authSeqId == residueNumber) {
                target = &c;
                targetIndex = i;
            }
        }
    }
    if (!target) {
        out.warnings.push_back("Chain " + chainId + " residue " + std::to_string(residueNumber) +
                               " (author numbering) is not in this structure.");
        return out;
    }
    out.wildType = target->residues[targetIndex].oneLetter();

    const std::string mutantName = threeLetterOf(mutant);
    if (mutantName == "UNK") {
        out.warnings.push_back(std::string("'") + mutant +
                               "' is not one of the 20 standard residues.");
        return out;
    }
    if (mutantName == "PRO") {
        out.warnings.push_back(
            "Proline is refused: its side chain closes a ring back onto the backbone nitrogen, "
            "so it constrains phi and cannot be placed by side-chain rotamer construction on a "
            "fixed backbone. A proline substitution needs backbone modelling, which this build "
            "does not do.");
        return out;
    }
    if (mutantName == "GLY") {
        out.assumptions.push_back(
            "Glycine has no side chain, so nothing was built and there is nothing to repack.");
        out.rotamerProbability = 1.0;
        return out;
    }

    const bio::ResidueRotamers* tmpl = lib->residue(mutantName);
    if (!tmpl || tmpl->build.empty()) {
        out.warnings.push_back("The rotamer pack has no build tree for " + mutantName + ".");
        return out;
    }

    // --- assemble the flexible sites: the mutated position plus its neighbours ---
    std::vector<Site> sites;
    const auto siteFor = [&](const bio::Chain& c, std::size_t idx, const std::string& resName) {
        Site st;
        st.chain = &c;
        st.index = idx;
        st.residueName = resName;
        st.label = c.id + ":" + std::to_string(c.residues[idx].authSeqId);
        double phi = 0, psi = 0;
        const bool haveBackbone = bio::backboneTorsions(c, idx, phi, psi);
        if (haveBackbone) {
            st.rotamers = lib->rotamersAt(resName, phi, psi, st.backboneIndependent);
        } else {
            const bio::ResidueRotamers* rr = lib->residue(resName);
            if (rr) st.rotamers = rr->backboneIndependent;
            st.backboneIndependent = true;
        }
        return st;
    };

    Site mutated = siteFor(*target, targetIndex, mutantName);
    if (mutated.rotamers.empty()) {
        // ALA has a build tree (CB) but no chi and hence no rotamer rows.
        bio::RotamerEntry only;
        only.name = "-";
        only.probability = 1.0;
        only.count = tmpl->build.front().observations;
        mutated.rotamers.push_back(only);
    }
    if (mutated.backboneIndependent)
        out.assumptions.push_back(
            "The phi/psi bin at this position had fewer than the pack's minimum observations "
            "(or the residue is terminal), so the backbone-INDEPENDENT rotamer set was used.");
    sites.push_back(std::move(mutated));

    const bio::Atom* anchor = anchorAtom(target->residues[targetIndex]);
    if (!anchor) {
        out.warnings.push_back("The target residue has neither CB nor CA, so no neighbourhood "
                               "could be defined.");
        return out;
    }
    struct Cand {
        double d2 = 0;
        const bio::Chain* chain = nullptr;
        std::size_t index = 0;
    };
    std::vector<Cand> cands;
    for (const auto& c : model->chains) {
        for (std::size_t i = 0; i < c.residues.size(); ++i) {
            if (&c == target && i == targetIndex) continue;
            const bio::Residue& r = c.residues[i];
            const std::string& n = r.name;
            if (n == "GLY" || n == "ALA" || n == "PRO") continue;   // nothing to repack
            if (!lib->residue(n) || lib->residue(n)->backboneIndependent.empty()) continue;
            const bio::Atom* a = anchorAtom(r);
            if (!a) continue;
            const double dx = a->x - anchor->x, dy = a->y - anchor->y, dz = a->z - anchor->z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= kRepackRadiusAngstrom * kRepackRadiusAngstrom)
                cands.push_back({d2, &c, i});
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.d2 < b.d2; });
    if (static_cast<int>(cands.size()) > kMaxRepackedNeighbours)
        cands.resize(static_cast<std::size_t>(kMaxRepackedNeighbours));
    for (const Cand& c : cands)
        sites.push_back(siteFor(*c.chain, c.index, c.chain->residues[c.index].name));

    // --- build every rotamer of every site, and collect the fixed environment ---
    for (Site& st : sites) {
        const bio::ResidueRotamers* rr = lib->residue(st.residueName);
        for (const bio::RotamerEntry& e : st.rotamers)
            st.built.push_back(rr ? bio::buildSideChain(*rr, st.chain->residues[st.index], e.chi)
                                  : std::vector<bio::Atom>{});
    }

    // Pointers into the structure, not copies: the exclusion test below is an
    // identity test, so a copied atom would never match its own exclusion entry.
    std::vector<const bio::Atom*> fixed;
    for (const auto& c : model->chains) {
        for (std::size_t i = 0; i < c.residues.size(); ++i) {
            bool flexible = false;
            for (const Site& st : sites)
                if (st.chain == &c && st.index == i) flexible = true;
            for (const bio::Atom& a : c.residues[i].atoms) {
                if (a.element == "H" || a.element == "D") continue;
                if (a.hetatm) continue;   // waters and ligands are not the packing environment
                if (flexible && !isBackboneName(a.name)) continue;   // its side chain is rebuilt
                fixed.push_back(&a);
            }
        }
    }
    // Backbone atoms of each flexible site and of its sequence neighbours are
    // 1-2/1-3/1-4 bonded to the side chain being placed; they sit inside the clash
    // cut-off for reasons no rotamer choice can change, so they are excluded.
    const auto excludedFor = [&](const Site& st) {
        std::vector<const bio::Atom*> ex;
        const auto add = [&ex](const bio::Residue& r) {
            for (const bio::Atom& a : r.atoms)
                if (isBackboneName(a.name)) ex.push_back(&a);
        };
        add(st.chain->residues[st.index]);
        if (st.index > 0) add(st.chain->residues[st.index - 1]);
        if (st.index + 1 < st.chain->residues.size()) add(st.chain->residues[st.index + 1]);
        return ex;
    };

    bio::DeeProblem problem;
    const std::size_t np = sites.size();
    problem.self.resize(np);
    problem.pair.assign(np, std::vector<std::vector<std::vector<double>>>(np));

    const auto environmentFor = [&](const Site& st) {
        const std::vector<const bio::Atom*> ex = excludedFor(st);
        std::vector<bio::Atom> env;
        env.reserve(fixed.size());
        for (const bio::Atom* a : fixed) {
            if (std::find(ex.begin(), ex.end(), a) != ex.end()) continue;
            env.push_back(*a);
        }
        return env;
    };

    for (std::size_t i = 0; i < np; ++i) {
        const std::vector<bio::Atom> env = environmentFor(sites[i]);
        problem.self[i].resize(sites[i].built.size());
        for (std::size_t r = 0; r < sites[i].built.size(); ++r) {
            const int clashes = bio::countClashes(sites[i].built[r], env);
            const double p = std::max(1e-6, sites[i].rotamers[r].probability);
            // Unit-free selection objective: one clash costs one, and an unlikely
            // rotamer costs -ln(p). Never reported; RotamerRebuild has no energy.
            problem.self[i][r] = static_cast<double>(clashes) - std::log(p);
        }
    }
    for (std::size_t i = 0; i < np; ++i) {
        for (std::size_t j = 0; j < np; ++j) {
            problem.pair[i][j].assign(sites[i].built.size(),
                                      std::vector<double>(sites[j].built.size(), 0.0));
            if (i == j) continue;
            for (std::size_t r = 0; r < sites[i].built.size(); ++r)
                for (std::size_t t = 0; t < sites[j].built.size(); ++t)
                    problem.pair[i][j][r][t] =
                        static_cast<double>(bio::countClashes(sites[i].built[r], sites[j].built[t]));
        }
    }

    const bio::DeeResult dee = bio::goldsteinDee(problem);
    const int pick = dee.chosen.empty() ? 0 : dee.chosen[0];
    const bio::RotamerEntry& chosen = sites[0].rotamers[static_cast<std::size_t>(pick)];
    out.chiAngles = chosen.chi;
    out.rotamerProbability = chosen.probability;

    // The reported clash count is the mutated side chain against everything else in
    // the final packing: the fixed environment plus the repacked neighbours.
    {
        const std::vector<bio::Atom> env = environmentFor(sites[0]);
        int clashes = bio::countClashes(sites[0].built[static_cast<std::size_t>(pick)], env);
        for (std::size_t j = 1; j < np; ++j)
            clashes += bio::countClashes(sites[0].built[static_cast<std::size_t>(pick)],
                                         sites[j].built[static_cast<std::size_t>(dee.chosen[j])]);
        out.clashCount = clashes;
    }
    for (std::size_t j = 1; j < np; ++j) out.repackedNeighbours.push_back(sites[j].label);

    out.assumptions.push_back(
        "The backbone is held fixed. Only the mutated side chain and " +
        std::to_string(np - 1) + " neighbouring side chain(s) within " +
        std::to_string(static_cast<int>(kRepackRadiusAngstrom)) +
        " A (CB-CB, nearest first, capped at " + std::to_string(kMaxRepackedNeighbours) +
        ") were allowed to move.");
    out.assumptions.push_back(
        "Rotamers were selected by Goldstein dead-end elimination over heavy-atom clash counts "
        "and library probability only, followed by " +
        std::string(dee.exhaustive ? "an exhaustive" : "an iterative single-position") +
        " repack of the survivors (" + std::to_string(dee.eliminated) +
        " rotamers eliminated in " + std::to_string(dee.passes) + " passes).");
    out.assumptions.push_back(
        "A clash is a heavy-atom pair closer than " +
        std::to_string(bio::kClashOverlapFactor).substr(0, 4) +
        " x the sum of their van der Waals radii. Hydrogens, waters and ligands are not in the "
        "environment, and bonded backbone atoms of the packed residue and its sequence "
        "neighbours are excluded.");
    out.assumptions.push_back(
        "This is a constructed side chain, not an energy calculation. There is no force field "
        "behind it and no stability, affinity or activity claim follows from it.");
    return out;
}

}  // namespace biocad
