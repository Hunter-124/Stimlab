// Functional-group perception: the SMARTS rule pack in
// assets/packs/rules/functional-groups.json, flag by flag.
//
// WHY the table is written negative-first: every flag here feeds an ADMET,
// stability or bioavailability verdict, so a FALSE POSITIVE is a wrong claim
// about a real drug, not a harmless extra. The cases that matter most are
// therefore the look-alikes - amphetamine's alpha-methyl versus phenethylamine's
// methylene, a benzylamine versus a phenethylamine, a ring amine versus a
// flexible one, a lactam versus an anilide, a quaternary ammonium versus a base.
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "chem/Analysis.h"
#include "chem/Smiles.h"

using namespace biocad::chem;

namespace {

FunctionalGroups groupsOf(const std::string& smiles) {
    auto m = parseSmiles(smiles);
    REQUIRE(m.has_value());
    return detectGroups(*m);
}

}  // namespace

// A rule that failed to parse leaves its flag permanently false, which is
// indistinguishable from the flag being honestly false - so the pack loading
// clean is itself an assertion, and it runs first.
TEST_CASE("The functional-group rule pack loads with no errors", "[chem][groups]") {
    for (const auto& e : groupPackErrors()) WARN(e);
    REQUIRE(groupPackErrors().empty());
}

TEST_CASE("Functional-group flags over the shipped library", "[chem][groups]") {
    struct Case {
        const char* id;
        const char* smiles;
        std::vector<bool FunctionalGroups::*> on;   // must be set
        std::vector<bool FunctionalGroups::*> off;  // must NOT be set
    };

    constexpr auto kArom = &FunctionalGroups::aromaticRing;
    constexpr auto kPhenol = &FunctionalGroups::phenol;
    constexpr auto kCatechol = &FunctionalGroups::catechol;
    constexpr auto kEster = &FunctionalGroups::ester;
    constexpr auto kAcid = &FunctionalGroups::carboxylicAcid;
    constexpr auto kAmide = &FunctionalGroups::amide;
    constexpr auto kKetone = &FunctionalGroups::ketone;
    constexpr auto kArylKetone = &FunctionalGroups::arylKetone;
    constexpr auto kAldehyde = &FunctionalGroups::aldehyde;
    constexpr auto kEther = &FunctionalGroups::ether;
    constexpr auto kMdo = &FunctionalGroups::methylenedioxy;
    constexpr auto kPrimary = &FunctionalGroups::primaryAmine;
    constexpr auto kSecondary = &FunctionalGroups::secondaryAmine;
    constexpr auto kTertiary = &FunctionalGroups::tertiaryAmine;
    constexpr auto kBasic = &FunctionalGroups::basicAmine;
    constexpr auto kNitrile = &FunctionalGroups::nitrile;
    constexpr auto kNitro = &FunctionalGroups::nitro;
    constexpr auto kHalogen = &FunctionalGroups::halogen;
    constexpr auto kSulfoxide = &FunctionalGroups::sulfoxide;
    constexpr auto kSulfone = &FunctionalGroups::sulfone;
    constexpr auto kPea = &FunctionalGroups::phenethylamine;
    constexpr auto kCatecholamine = &FunctionalGroups::catecholamine;
    constexpr auto kAnilide = &FunctionalGroups::anilide;
    constexpr auto kMao = &FunctionalGroups::maoLabileAmine;

    const std::vector<Case> cases = {
        // -- the MAO discrimination this pack exists to get right --------------
        // 2-phenylethylamine: unsubstituted alpha methylene, MAO-B's canonical
        // substrate. Oral bioavailability is essentially nil because of it.
        {"phenethylamine", "NCCc1ccccc1", {kArom, kPrimary, kBasic, kPea, kMao}, {kPhenol, kAnilide}},
        // Amphetamine: the alpha carbon carries a methyl, so it is a CH and MAO
        // cannot deaminate it. This single H is the difference between a drug and
        // a first-pass casualty, and it must not be flagged.
        {"amphetamine", "CC(N)Cc1ccccc1", {kArom, kPrimary, kBasic, kPea},
         {kMao, kEster, kAnilide, kCatechol}},
        // Methamphetamine: N-methyl (a CH3, three hydrogens, not the alpha
        // methylene) AND an alpha methyl. Still not MAO-labile.
        {"methamphetamine", "CNC(C)Cc1ccccc1", {kArom, kSecondary, kBasic, kPea},
         {kMao, kPrimary, kTertiary}},
        {"tyramine", "NCCc1ccc(O)cc1", {kArom, kPhenol, kPrimary, kBasic, kPea, kMao},
         {kCatechol, kCatecholamine}},
        {"dopamine", "NCCc1ccc(O)c(O)c1",
         {kArom, kPhenol, kCatechol, kPrimary, kBasic, kPea, kCatecholamine, kMao}, {kEther}},
        {"norepinephrine", "NCC(O)c1ccc(O)c(O)c1",
         {kPhenol, kCatechol, kPrimary, kPea, kCatecholamine, kMao}, {kEther, kAnilide}},
        {"serotonin", "NCCc1c[nH]c2ccc(O)cc12", {kArom, kPhenol, kPrimary, kBasic, kPea, kMao},
         {kCatechol, kAmide}},
        // Nicotine's amine is a ring (pyrrolidine) nitrogen one bond from the
        // arene: neither MAO-labile nor a phenethylamine core.
        {"nicotine", "CN1CCCC1c1cccnc1", {kArom, kTertiary, kBasic}, {kMao, kPea, kPrimary}},
        // Cocaine: tropane ring amine plus two esters. Not MAO-labile, not a
        // phenethylamine.
        {"cocaine", "COC(=O)C1C(OC(=O)c2ccccc2)CC2CCC1N2C", {kArom, kEster, kTertiary, kBasic},
         {kMao, kPea, kAcid, kPrimary}},
        // MDMA: alpha-methyl again, so not MAO-labile, but the methylenedioxy
        // bridge (CYP2D6 demethylenation) must be seen.
        {"mdma", "CNC(C)Cc1ccc2OCOc2c1", {kArom, kEther, kMdo, kSecondary, kBasic, kPea},
         {kMao, kPhenol, kCatechol}},
        {"methylphenidate", "COC(=O)C(c1ccccc1)C1CCCCN1", {kArom, kEster, kSecondary, kBasic, kPea},
         {kMao, kAcid}},

        // -- anilide versus lactam versus plain amide --------------------------
        // Acetaminophen: the para-aminophenol/anilide that bioactivates to NAPQI.
        {"acetaminophen", "CC(=O)Nc1ccc(O)cc1", {kArom, kPhenol, kAmide, kAnilide},
         {kBasic, kPrimary, kEster, kMao}},
        // Caffeine, written Kekule exactly as the pack authors it: the xanthine
        // carbonyls are amides and none of its nitrogens is a base (conjugate
        // acid pKa ~0.6). An N-acyl AROMATIC nitrogen is a lactam, not an
        // anilide, and does not form a quinone-imine.
        {"caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", {kArom, kAmide},
         {kBasic, kTertiary, kAnilide, kPhenol, kKetone}},
        {"theophylline", "Cn1c(=O)c2[nH]cnc2n(C)c1=O", {kArom, kAmide}, {kBasic, kAnilide, kKetone}},
        // Indomethacin's nitrogen is an N-acyl indole - aromatic, so again a
        // lactam-type amide rather than an anilide.
        {"indomethacin", "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
         {kArom, kAmide, kAcid, kEther, kHalogen}, {kAnilide, kBasic, kEster}},

        // -- acids, esters and the carboxylate that is not an ester ------------
        {"ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", {kArom, kAcid}, {kEster, kBasic, kPhenol}},
        {"aspirin", "CC(=O)Oc1ccccc1C(=O)O", {kArom, kAcid, kEster}, {kBasic, kPhenol, kKetone}},
        // L-carnitine is a zwitterion: a carboxylATE (an acid, not an ester) and
        // a quaternary ammonium (permanently charged, so not a basic amine at
        // all - it has no lone pair to protonate).
        {"l-carnitine", "C[N+](C)(C)CC(O)CC(=O)[O-]", {kAcid}, {kEster, kBasic, kTertiary, kArom}},

        // -- nitrogen basicity ------------------------------------------------
        // Aniline: a primary amine by structure, but its lone pair is conjugated
        // into the ring (conjugate acid pKa 4.6), so it is not a basic amine.
        {"aniline", "Nc1ccccc1", {kArom, kPrimary}, {kBasic, kSecondary, kMao}},
        // Diclofenac's nitrogen is a diarylamine - even less basic than aniline.
        {"diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", {kArom, kAcid, kSecondary, kHalogen},
         {kBasic, kMao, kPea}},
        {"diethylamine", "CCNCC", {kSecondary, kBasic}, {kArom, kPrimary, kTertiary}},
        {"triethylamine", "CCN(CC)CC", {kTertiary, kBasic}, {kArom, kPrimary, kSecondary}},
        // A sulfonamide N-H is acidic, not basic; celecoxib is not a base.
        {"celecoxib", "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F",
         {kArom, kSulfone, kHalogen, kPrimary}, {kBasic, kAcid}},
        // Cationic AROMATIC nitrogen: NMN's pyridinium and berberine's
        // isoquinolinium carry three sigma bonds and no hydrogen, which is exactly
        // the shape of a tertiary amine - but they are aromatic ring nitrogens
        // with no lone pair left, so neither is an amine of any class. These two
        // are the reason chem::perceiveAromaticity scores an N(+) in a pi ring
        // bond as a one-electron pyridine-type centre; scored as pyrrole-type the
        // rings failed Huckel, went aliphatic, and every N-methylated cofactor
        // read as a basic tertiary amine.
        {"nmn", "NC(=O)c1ccc[n+](c1)C1OC(COP(=O)(O)[O-])C(O)C1O", {kArom, kAmide, kEther},
         {kBasic, kTertiary, kPrimary, kAnilide}},
        {"berberine", "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2OCOc2cc1-3", {kArom, kEther, kMdo},
         {kBasic, kTertiary, kPea, kAmide}},

        // -- carbonyl classes, sulfur oxidation state, and the small groups ----
        {"benzaldehyde", "O=Cc1ccccc1", {kArom, kAldehyde}, {kKetone, kArylKetone, kAcid, kEster}},
        {"acetophenone", "CC(=O)c1ccccc1", {kArom, kKetone, kArylKetone}, {kAldehyde, kAcid}},
        {"cathinone", "CC(N)C(=O)c1ccccc1", {kArom, kKetone, kArylKetone, kPrimary, kBasic},
         {kAldehyde, kMao}},
        {"modafinil", "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1", {kArom, kAmide, kSulfoxide},
         {kSulfone, kBasic, kAnilide}},
        {"methanesulfonamide", "CS(=O)(=O)N", {kSulfone}, {kSulfoxide, kBasic, kArom}},
        {"benzonitrile", "N#Cc1ccccc1", {kArom, kNitrile}, {kBasic, kPrimary, kAmide}},
        {"nitrobenzene", "O=[N+]([O-])c1ccccc1", {kArom, kNitro},
         {kBasic, kTertiary, kPrimary, kNitrile}},
        {"anisole", "COc1ccccc1", {kArom, kEther}, {kMdo, kPhenol, kEster}},
        {"fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
         {kArom, kEther, kSecondary, kBasic, kHalogen}, {kMao, kPea, kPhenol}},
        // Codeine is written Kekule in the pack; the morphinan skeleton really
        // does contain an Ar-C-C-N core, which is why it is a CYP2D6 substrate.
        {"codeine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O",
         {kArom, kEther, kTertiary, kBasic, kPea}, {kMao, kPhenol, kEster}},
    };

    for (const auto& c : cases) {
        const auto g = groupsOf(c.smiles);
        for (const auto flag : c.on) {
            INFO(c.id << " should carry a flag it does not");
            REQUIRE(g.*flag);
        }
        for (const auto flag : c.off) {
            INFO(c.id << " carries a flag it must not");
            REQUIRE_FALSE(g.*flag);
        }
    }
    REQUIRE(cases.size() >= 20);
}

// Perception must not depend on how the SMILES happened to be spelled: the same
// molecule as a Kekule structure and as lowercase-aromatic must give the same
// flags. This is what the aromaticity pass buys, and it is why the pack's
// carbonyl patterns are written [#6X3] rather than [CX3].
TEST_CASE("Group perception is independent of Kekule vs aromatic spelling", "[chem][groups]") {
    struct Pair { const char* kekule; const char* aromatic; };
    const Pair pairs[] = {
        {"CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"},  // caffeine
        {"CC(N)CC1=CC=CC=C1", "CC(N)Cc1ccccc1"},                          // amphetamine
        {"NCCC1=CC(O)=C(O)C=C1", "NCCc1cc(O)c(O)cc1"},                    // dopamine
    };
    for (const auto& p : pairs) {
        const auto a = groupsOf(p.kekule);
        const auto b = groupsOf(p.aromatic);
        INFO(p.kekule << " vs " << p.aromatic);
        REQUIRE(a.aromaticRing == b.aromaticRing);
        REQUIRE(a.amide == b.amide);
        REQUIRE(a.basicAmine == b.basicAmine);
        REQUIRE(a.tertiaryAmine == b.tertiaryAmine);
        REQUIRE(a.phenethylamine == b.phenethylamine);
        REQUIRE(a.catechol == b.catechol);
        REQUIRE(a.catecholamine == b.catecholamine);
        REQUIRE(a.maoLabileAmine == b.maoLabileAmine);
    }
}

TEST_CASE("An empty molecule carries no flags", "[chem][groups]") {
    const Molecule empty;
    const auto g = detectGroups(empty);
    REQUIRE_FALSE(g.aromaticRing);
    REQUIRE_FALSE(g.basicAmine);
    REQUIRE_FALSE(g.anilide);
}
