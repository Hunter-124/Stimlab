#include "modules/BiologicsModule.h"

#include <algorithm>

#include "bio/Imgt.h"
#include "bio/Interface.h"
#include "bio/Liabilities.h"

namespace biocad {
namespace {

// Which chain of a structure a numbered domain came from, matched by sequence so a
// caller does not have to pass a chain id it may not know.
std::string chainCarrying(const AbDomain& d, const bio::Structure* s) {
    if (!s || d.sequence.size() < 20) return {};
    const bio::Model* m = s->model(1);
    if (!m) return {};
    const std::string probe = d.sequence.substr(0, 20);
    for (const auto& c : m->chains) {
        const std::vector<char> sv = bio::sequenceOf(c);
        if (std::string(sv.begin(), sv.end()).find(probe) != std::string::npos) return c.id;
    }
    return {};
}

}  // namespace

AbDomain RealBiologics::number(const std::string& sequence, NumberingScheme scheme) const {
    // IMGT is canonical: everything is numbered in IMGT first and then VIEWED in the
    // requested scheme, so a scheme with no published table cannot produce numbering
    // by accident.
    const AbDomain imgt = bio::numberDomain(sequence);
    if (scheme == NumberingScheme::Imgt || !imgt.numbered) return imgt;
    return bio::convertScheme(imgt, scheme);
}

AbDomain RealBiologics::convertScheme(const AbDomain& domain, NumberingScheme to) const {
    return bio::convertScheme(domain, to);
}

std::vector<SequenceLiability> RealBiologics::liabilities(const AbDomain& domain,
                                                         const bio::Structure* structure) const {
    return bio::scanLiabilities(domain, structure, chainCarrying(domain, structure));
}

DevelopabilityReport RealBiologics::developability(const std::vector<std::string>& chains,
                                                  const bio::Structure* structure) const {
    bio::DevelopabilityInput in;
    in.chains = chains;
    in.structure = structure;
    if (structure) {
        // The origin is what the structure says it is. It is NOT invented here: an
        // empty entry id leaves it empty, which keeps the TAP metrics NotComputed.
        if (!structure->id.empty()) in.structureOrigin = "structure " + structure->id + ", as " +
                                                         "loaded from " + structure->source;
        for (const auto& c : chains) {
            AbDomain d = bio::numberDomain(c);
            d.sequence = c;
            in.chainIds.push_back(chainCarrying(d, structure));
        }
    }
    return bio::developability(in);
}

MassLadder RealBiologics::massLadder(const std::vector<std::string>& chains,
                                     int disulfideCount) const {
    bio::MassLadderInput in;
    for (const auto& c : chains) {
        const AbDomain d = bio::numberDomain(c);
        const bool heavy = d.numbered && (d.chain == AbChainType::HeavyVh ||
                                          d.chain == AbChainType::Vhh);
        (heavy ? in.heavyChains : in.lightChains).push_back(c);
    }
    if (in.heavyChains.empty() && in.lightChains.size() == chains.size()) {
        // Nothing numbered as a heavy chain: treat the list as given rather than
        // inventing a glycosylation site on a chain that may not have one.
        in.heavyChains = in.lightChains;
        in.lightChains.clear();
        in.includeGlycoforms = false;
    }
    in.intrachainDisulfides = std::max(0, disulfideCount);
    MassLadder out = bio::massLadder(in);
    out.assumptions.push_back(
        "Heavy/light assignment: " + std::to_string(in.heavyChains.size()) +
        " chain(s) numbered as VH or VHH and are treated as heavy, " +
        std::to_string(in.lightChains.size()) +
        " as light. The glycoform ladder puts one glycan on each heavy chain; if the "
        "construct is aglycosylated, read the intact row only." +
        (in.includeGlycoforms ? "" : " No glycoform row was produced."));
    return out;
}

PeptideMap RealBiologics::digest(const std::string& chain, const std::string& protease,
                                 int maxMissedCleavages) const {
    bio::DigestOptions o;
    o.protease = protease;
    o.maxMissedCleavages = std::max(0, maxMissedCleavages);
    // The two modifications a biologics peptide map is actually read for. Both are
    // formula arithmetic, so their deltas come from the NIST table like every mass.
    o.mods = {{"oxidation", "O", "", "MW"}, {"deamidation", "O", "NH", "N"}};
    return bio::digest(chain, o);
}

InterfaceReport RealBiologics::interfaceOf(const bio::Structure& complex,
                                           const std::string& chainsA,
                                           const std::string& chainsB) const {
    bio::InterfaceOptions o;
    // Which chains are the antibody is DERIVED, by numbering them: a chain that
    // numbers as a V-DOMAIN gets its CDR contacts reported. A chain that does not is
    // simply the antigen side, never a guess from its letter.
    const bio::Model* m = complex.model(1);
    if (m) {
        const auto a = bio::parseChainList(chainsA), b = bio::parseChainList(chainsB);
        for (const auto& c : m->chains) {
            const bool named = std::find(a.begin(), a.end(), c.id) != a.end() ||
                               std::find(b.begin(), b.end(), c.id) != b.end();
            if (!named) continue;
            const std::vector<char> sv = bio::sequenceOf(c);
            if (bio::numberDomain(std::string(sv.begin(), sv.end())).numbered)
                o.antibodyChains.push_back(c.id);
        }
    }
    return bio::interfaceOf(complex, chainsA, chainsB, o);
}

AlanineScanReport RealBiologics::alanineScan(const bio::Structure& complex,
                                             const std::string& chainsA,
                                             const std::string& chainsB) const {
    return bio::alanineScan(complex, chainsA, chainsB);
}

}  // namespace biocad
