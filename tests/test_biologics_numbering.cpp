// IMGT numbering, chain typing, insertion codes, refusals and scheme conversion.
//
// The sequences are REAL and named: they were taken from the RCSB entry FASTA of the
// listed PDB entries (https://www.rcsb.org/fasta/entry/<id>), so anyone can re-fetch
// them and check the CDRs by eye.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>

#include "bio/Imgt.h"

using namespace biocad;
using namespace biocad::bio;

namespace {

// 1N8Z: trastuzumab (Herceptin) Fab bound to the HER2 ectodomain. Chain B is the
// heavy chain (VH + CH1), chain A the complete light chain (VL + CL).
constexpr const char* kTrastuzumabH =
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAED"
    "TAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSS"
    "GLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEP";
constexpr const char* kTrastuzumabL =
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQ"
    "HYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLT"
    "LSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC";
// 1VFB: anti-lysozyme D1.3 Fv, chain A (kappa light).
constexpr const char* kD13L =
    "DIVLTQSPASLSASVGETVTITCRASGNIHNYLAWYQQKQGKSPQLLVYYTTTLADGVPSRFSGSGSGTQYSLKINSLQPEDFGSYYCQH"
    "FWSTPRTFGGGTKLEIK";
// 2FB4: IgG1-lambda Kol Fab, chain A[auth L] (lambda light).
constexpr const char* kKolL =
    "QSVLTQPPSASGTPGQRVTISCSGTSSNIGSSTVNWYQQLPGMAPKLLIYRDAMRPSGVPDRFSGSKSGASASLAIGGLQSEDETDYYCA"
    "AWDVSLNAYVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASS"
    "YLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS";
// 1MEL: camelid VHH against lysozyme, chain A. Its CDR3 is long - the insertion case.
constexpr const char* kVhh =
    "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPED"
    "TAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSSGRYPYDVPDYGSGRA";
// 1TCR: mouse alpha/beta T-cell receptor, chain B (beta).
constexpr const char* kTcrBeta =
    "EAAVTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYF"
    "CASGGGGTLYFGAGTRLSVLEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSY"
    "CLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADC";

char residueAt(const AbDomain& d, int position) {
    for (const auto& r : d.residues)
        if (r.position == position && r.insertionCode.empty()) return r.aminoAcid;
    return '?';
}

std::string regionSequence(const AbDomain& d, const std::string& region) {
    std::string s;
    for (const auto& r : d.residues)
        if (r.region == region) s.push_back(r.aminoAcid);
    return s;
}

void requireAllFiveAnchors(const AbDomain& d) {
    REQUIRE(d.numbered);
    REQUIRE(d.anchorFailures.empty());
    CHECK(residueAt(d, 23) == 'C');    // 1st-CYS
    CHECK(residueAt(d, 41) == 'W');    // CONSERVED-TRP
    CHECK(std::string("AVLIMFCWYP").find(residueAt(d, 89)) != std::string::npos);
    CHECK(residueAt(d, 104) == 'C');   // 2nd-CYS
    const char j = residueAt(d, 118);
    CHECK((j == 'F' || j == 'W'));     // J-PHE or J-TRP
}

}  // namespace

TEST_CASE("the IMGT reference pack loads", "[biologics]") {
    const ImgtReference& ref = imgtReference();
    for (const auto& e : ref.errors) WARN(e);
    REQUIRE(ref.ok);
    CHECK(ref.vGenes.size() > 300);
    CHECK(ref.jGenes.size() > 20);
    CHECK(ref.licence.find("CC BY-4.0") != std::string::npos);
    // Index 104 of every gapped V profile is the 2nd-CYS. That invariant is what makes
    // position transfer exact rather than inferred, so it is asserted over the pack.
    for (const auto& g : ref.vGenes) {
        REQUIRE(g.gapped.size() >= 104);
        REQUIRE(g.gapped[103] == 'C');
    }
}

TEST_CASE("a published VH numbers with all five anchors and its known CDR lengths",
          "[biologics]") {
    const AbDomain d = numberDomain(kTrastuzumabH);
    requireAllFiveAnchors(d);
    CHECK(d.chain == AbChainType::HeavyVh);
    // Trastuzumab's IMGT CDRs, readable straight off the sequence above.
    CHECK(regionSequence(d, "CDR1") == "GFNIKDTY");
    CHECK(regionSequence(d, "CDR2") == "IYPTNGYT");
    CHECK(regionSequence(d, "CDR3") == "SRWGGDGFYAMDY");
    REQUIRE(d.cdrLengths.size() == 3);
    CHECK(d.cdrLengths[0] == 8);
    CHECK(d.cdrLengths[1] == 8);
    CHECK(d.cdrLengths[2] == 13);
    CHECK(regionSequence(d, "FR4") == "WGQGTLVTVSS");
    // The chain carries CH1: everything past IMGT 128 is refused with a warning
    // rather than numbered as if the constant domain were part of the V-DOMAIN.
    CHECK_FALSE(d.warnings.empty());
    CHECK(d.closestGermlineSet.rfind("IGHV", 0) == 0);
}

TEST_CASE("kappa and lambda light chains number and are typed apart", "[biologics]") {
    const AbDomain vk = numberDomain(kTrastuzumabL);
    requireAllFiveAnchors(vk);
    CHECK(vk.chain == AbChainType::LightVKappa);
    CHECK(regionSequence(vk, "CDR1") == "QDVNTA");
    CHECK(regionSequence(vk, "CDR2") == "SAS");
    CHECK(regionSequence(vk, "CDR3") == "QQHYTTPPT");

    const AbDomain d13 = numberDomain(kD13L);
    requireAllFiveAnchors(d13);
    CHECK(d13.chain == AbChainType::LightVKappa);
    CHECK(regionSequence(d13, "CDR3") == "QHFWSTPRT");

    const AbDomain lambda = numberDomain(kKolL);
    requireAllFiveAnchors(lambda);
    CHECK(lambda.chain == AbChainType::LightVLambda);
    CHECK(lambda.closestGermlineSet.rfind("IGLV", 0) == 0);
}

TEST_CASE("a long CDR3 numbers with insertion codes and its neighbours keep their positions",
          "[biologics]") {
    const AbDomain d = numberDomain(kVhh);
    requireAllFiveAnchors(d);
    CHECK(d.chain == AbChainType::Vhh);
    REQUIRE(d.cdrLengths.size() == 3);
    CHECK(d.cdrLengths[2] == 26);   // longer than the 13 IMGT CDR3 slots

    int insertions = 0;
    for (const auto& r : d.residues)
        if (!r.insertionCode.empty()) ++insertions;
    CHECK(insertions == d.cdrLengths[2] - 13);

    // Insertions land on 111 and 112 ONLY, which is what keeps every other CDR3
    // position where a 13-residue loop would have put it.
    for (const auto& r : d.residues)
        if (!r.insertionCode.empty()) CHECK((r.position == 111 || r.position == 112));
    for (int p : {105, 106, 107, 108, 109, 110, 113, 114, 115, 116, 117})
        CHECK(residueAt(d, p) != '?');

    // The neighbours of the insertion point are unshifted: position order across the
    // whole domain is still monotonic in sequence index.
    int lastPos = 0;
    for (const auto& r : d.residues) {
        CHECK(r.position >= lastPos);
        lastPos = r.position;
    }
}

TEST_CASE("a T-cell receptor beta chain is rejected with its reason listed", "[biologics]") {
    const AbDomain d = numberDomain(kTcrBeta);
    CHECK(d.chain == AbChainType::TcrBeta);
    CHECK_FALSE(d.numbered);
    CHECK(d.residues.empty());          // no numbering at all, not a partial one
    REQUIRE_FALSE(d.anchorFailures.empty());
    CHECK(d.anchorFailures.front().find("T-cell receptor") != std::string::npos);
    CHECK(d.cdrLengths.empty());
}

TEST_CASE("something far too short to be a V-DOMAIN is refused", "[biologics]") {
    const AbDomain d = numberDomain("EVQLVESGGGLVQPGGSLRLSCAAS");
    CHECK_FALSE(d.numbered);
    CHECK(d.residues.empty());
    REQUIRE_FALSE(d.anchorFailures.empty());
}

TEST_CASE("IMGT -> Kabat -> IMGT round-trips identically", "[biologics]") {
    for (const char* seq : {kTrastuzumabH, kTrastuzumabL, kD13L, kKolL, kVhh}) {
        const AbDomain imgt = numberDomain(seq);
        REQUIRE(imgt.numbered);
        const AbDomain kabat = convertScheme(imgt, NumberingScheme::Kabat);
        REQUIRE(kabat.numbered);
        REQUIRE(kabat.residues.size() == imgt.residues.size());
        CHECK(kabat.scheme == NumberingScheme::Kabat);
        const AbDomain back = convertScheme(kabat, NumberingScheme::Imgt);
        REQUIRE(back.numbered);
        REQUIRE(back.residues.size() == imgt.residues.size());
        for (std::size_t i = 0; i < imgt.residues.size(); ++i) {
            CHECK(back.residues[i].position == imgt.residues[i].position);
            CHECK(back.residues[i].insertionCode == imgt.residues[i].insertionCode);
            CHECK(back.residues[i].aminoAcid == imgt.residues[i].aminoAcid);
            CHECK(back.residues[i].region == imgt.residues[i].region);
        }
    }
}

TEST_CASE("the Kabat landmarks match the published correspondence", "[biologics]") {
    auto kabatOf = [](const AbDomain& imgt, int position) {
        int index = -1;
        for (const auto& r : imgt.residues)
            if (r.position == position && r.insertionCode.empty()) index = r.sequenceIndex;
        const AbDomain kabat = convertScheme(imgt, NumberingScheme::Kabat);
        for (const auto& r : kabat.residues)
            if (r.sequenceIndex == index) return r.position;
        return -1;
    };
    const AbDomain heavy = numberDomain(kTrastuzumabH);
    const AbDomain light = numberDomain(kTrastuzumabL);
    REQUIRE(heavy.numbered);
    REQUIRE(light.numbered);
    // IMGT 23/41/104 are Kabat H22/H36/H92 and L23/L35/L88.
    CHECK(kabatOf(heavy, 23) == 22);
    CHECK(kabatOf(heavy, 41) == 36);
    CHECK(kabatOf(heavy, 104) == 92);
    CHECK(kabatOf(light, 23) == 23);
    CHECK(kabatOf(light, 41) == 35);
    CHECK(kabatOf(light, 104) == 88);
}

TEST_CASE("a scheme with no obtainable table is refused, not guessed", "[biologics]") {
    const AbDomain imgt = numberDomain(kTrastuzumabH);
    REQUIRE(imgt.numbered);
    for (auto scheme : {NumberingScheme::Martin, NumberingScheme::Aho}) {
        const AbDomain out = convertScheme(imgt, scheme);
        CHECK_FALSE(out.numbered);
        CHECK(out.residues.empty());
        REQUIRE_FALSE(out.warnings.empty());
        CHECK(out.warnings.back().find("REFUSED") != std::string::npos);
    }
}
