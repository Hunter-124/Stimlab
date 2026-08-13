// DNA/RNA sequence core: the IUPAC alphabet, genetic codes, translation, ORFs,
// GenBank/FASTA I/O and restriction mapping.
//
// The load-bearing tests are the ones that catch arithmetic rather than plumbing:
// reverse complement as an involution over the whole ambiguity alphabet, the
// committed pUC19 record's CDSs translating to their own /translation qualifiers
// (which exercises complement() locations, multi-line quoted qualifiers and the
// genetic code at once), and the circular cut arithmetic - a site spanning the
// origin must be found exactly once, and fragment lengths must sum to the
// template length exactly.

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>

#include "bio/Codon.h"
#include "bio/NucIo.h"
#include "bio/NucSeq.h"
#include "bio/Restriction.h"
#include "core/Error.h"

using namespace biocad;
using namespace biocad::bio;

namespace {

std::filesystem::path fixture(const char* name) {
    return std::filesystem::path(BIOCAD_TEST_FIXTURES) / name;
}

std::filesystem::path pack(const char* name) {
    return std::filesystem::path(BIOCAD_ASSETS_DIR) / "packs" / "nucleic" / name;
}

std::string readAll(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

const NucRecord& puc19() {
    static const std::vector<NucRecord> records = readGenBank(readAll(fixture("pUC19.gb")));
    REQUIRE(records.size() == 1);
    return records[0];
}

}  // namespace

TEST_CASE("reverse complement is an involution over arbitrary IUPAC strings", "[bio][nuc]") {
    const std::string alphabet = "ACGTRYSWKMBDHVN";
    std::mt19937 rng(20260813);
    std::uniform_int_distribution<int> lenDist(1, 60);
    std::uniform_int_distribution<int> symDist(0, static_cast<int>(alphabet.size()) - 1);
    for (int i = 0; i < 1000; ++i) {
        std::string s;
        const int n = lenDist(rng);
        for (int k = 0; k < n; ++k) s.push_back(alphabet[static_cast<std::size_t>(symDist(rng))]);
        REQUIRE(reverseComplement(reverseComplement(s)) == s);
    }
    // The ambiguity pairings, spelled out, because this is where the mapping is
    // usually wrong: R<->Y, K<->M, B<->V, D<->H and S/W/N self-complementary.
    REQUIRE(reverseComplement("RYKMBVDHSWN") == "NWSDHBVKMRY");
    REQUIRE(complementBase('R') == 'Y');
    REQUIRE(complementBase('B') == 'V');
    REQUIRE(complementBase('D') == 'H');
    REQUIRE(complementBase('S') == 'S');
    REQUIRE(complementBase('W') == 'W');
    REQUIRE(complementBase('N') == 'N');
    REQUIRE(complementBase('-') == 'N');   // a gap is not silently promoted to a base
    REQUIRE(reverseComplementRna("ACGT") == "ACGU");
}

TEST_CASE("IUPAC matching is asymmetric on purpose", "[bio][nuc]") {
    REQUIRE(iupacMatches('N', 'A'));
    REQUIRE(iupacMatches('R', 'G'));
    REQUIRE_FALSE(iupacMatches('R', 'C'));
    // A pattern base must not match an ambiguous observation: a site is present
    // only where the sequence is known to contain it.
    REQUIRE_FALSE(iupacMatches('A', 'N'));
    REQUIRE_FALSE(iupacMatches('A', 'R'));
}

TEST_CASE("NCBI genetic code tables differ where the published tables differ", "[bio][nuc]") {
    const GeneticCodeTables tables = loadGeneticCodes(pack("genetic-codes.json"));
    for (int id : {1, 2, 3, 4, 5, 11}) REQUIRE(tables.find(id) != nullptr);
    const GeneticCode& t1 = *tables.find(1);
    const GeneticCode& t2 = *tables.find(2);
    const GeneticCode& t11 = *tables.find(11);
    REQUIRE(t1.translateCodon("TGA") == '*');
    REQUIRE(t2.translateCodon("TGA") == 'W');
    REQUIRE(t1.translateCodon("AGA") == 'R');
    REQUIRE(t2.translateCodon("AGA") == '*');
    REQUIRE(t1.translateCodon("ATA") == 'I');
    REQUIRE(t2.translateCodon("ATA") == 'M');
    REQUIRE(tables.find(3)->translateCodon("CTG") == 'T');
    REQUIRE(tables.find(5)->translateCodon("AGA") == 'S');
    REQUIRE(tables.find(4)->translateCodon("TGA") == 'W');
    // Table 11 assigns the same amino acids as table 1 and differs only in its
    // documented initiators.
    for (const char* c : {"TTT", "CTG", "ATA", "AGA", "TGA", "GGG"}) {
        REQUIRE(t11.translateCodon(c) == t1.translateCodon(c));
    }
    REQUIRE(t11.isStart("ATT"));
    REQUIRE_FALSE(t1.isStart("ATT"));
    REQUIRE(t1.isStart("ATG"));
    REQUIRE(t1.isStart("TTG"));
    // An ambiguous codon translates when every expansion agrees, and only then.
    REQUIRE(t1.translateCodon("GCN") == 'A');
    REQUIRE(t1.translateCodon("CTN") == 'L');
    REQUIRE(t1.translateCodon("AAN") == 'X');
}

TEST_CASE("the pUC19 fixture's CDSs translate to their own /translation", "[bio][nuc]") {
    const NucRecord& rec = puc19();
    REQUIRE(rec.sequence.size() == 2686);
    REQUIRE(rec.circular);
    REQUIRE(rec.features.size() == 4);
    const GeneticCodeTables tables = loadGeneticCodes(pack("genetic-codes.json"));

    int cdsCount = 0;
    for (const auto& feature : rec.features) {
        if (feature.type != "CDS") continue;
        ++cdsCount;
        std::string declared;
        for (const auto& q : feature.qualifiers) {
            if (q.first == "translation") declared = q.second;
        }
        REQUIRE_FALSE(declared.empty());
        // The multi-line /translation must have been rejoined without inserting a
        // space, which is the one qualifier where that matters.
        REQUIRE(declared.find(' ') == std::string::npos);
        REQUIRE(feature.strand == Strand::Reverse);
        std::string protein = translate(featureSequence(rec, feature), *tables.find(1));
        REQUIRE(protein.back() == '*');
        protein.pop_back();
        REQUIRE(protein == declared);
    }
    REQUIRE(cdsCount == 2);
    // TEM-1 beta-lactamase, as an independent anchor on the reading frame.
    for (const auto& feature : rec.features) {
        bool isBla = false;
        for (const auto& q : feature.qualifiers) {
            if (q.first == "gene" && q.second == "bla") isBla = true;
        }
        if (!isBla) continue;
        std::string protein = translate(featureSequence(rec, feature), *tables.find(1));
        protein.pop_back();
        REQUIRE(protein.size() == 286);
        REQUIRE(protein.compare(0, 23, "MSIQHFRVALIPFFAAFCLPVFA") == 0);
    }
}

TEST_CASE("GenBank locations cover the whole grammar", "[bio][nuc]") {
    std::vector<std::string> warnings;
    const ParsedLocation joined = parseLocation("complement(join(<12..20,30..40>))", &warnings);
    REQUIRE(joined.strand == Strand::Reverse);
    REQUIRE(joined.parts.size() == 2);
    REQUIRE(joined.parts[0] == std::make_pair(11, 20));
    REQUIRE(joined.parts[1] == std::make_pair(29, 40));
    REQUIRE(joined.partial5);
    REQUIRE(joined.partial3);
    REQUIRE_FALSE(warnings.empty());

    // The two spellings of a reverse-strand spliced feature must land on one
    // representation, or downstream code has to know which producer wrote the file.
    const ParsedLocation other = parseLocation("join(complement(30..40),complement(12..20))");
    REQUIRE(other.strand == Strand::Reverse);
    REQUIRE(other.parts.size() == 2);
    REQUIRE(other.parts[0] == std::make_pair(11, 20));
    REQUIRE(other.parts[1] == std::make_pair(29, 40));

    const ParsedLocation ordered = parseLocation("order(1..3,7..9)");
    REQUIRE(ordered.ordered);
    REQUIRE(ordered.parts.size() == 2);

    const ParsedLocation site = parseLocation("1023^1024");
    REQUIRE(site.betweenBases);
    REQUIRE(site.parts[0] == std::make_pair(1023, 1023));

    const ParsedLocation single = parseLocation("complement(7)");
    REQUIRE(single.parts[0] == std::make_pair(6, 7));
    REQUIRE(single.strand == Strand::Reverse);

    std::vector<std::string> remoteWarnings;
    const ParsedLocation remote = parseLocation("J00194.1:1..10", &remoteWarnings);
    REQUIRE(remote.parts.empty());
    REQUIRE_FALSE(remoteWarnings.empty());

    REQUIRE_THROWS_AS(parseLocation("not-a-location"), Error);

    REQUIRE(formatLocation({{11, 20}, {29, 40}}, Strand::Reverse) ==
            "complement(join(12..20,30..40))");
    REQUIRE(formatLocation({{0, 2686}}, Strand::Forward) == "1..2686");
}

TEST_CASE("a parsed GenBank record round-trips structurally", "[bio][nuc]") {
    const NucRecord& rec = puc19();
    const std::vector<NucRecord> again = readGenBank(writeGenBank(rec));
    REQUIRE(again.size() == 1);
    REQUIRE(again[0].id == rec.id);
    REQUIRE(again[0].description == rec.description);
    REQUIRE(again[0].circular == rec.circular);
    REQUIRE(again[0].sequence == rec.sequence);
    REQUIRE(again[0].features.size() == rec.features.size());
    for (std::size_t i = 0; i < rec.features.size(); ++i) {
        REQUIRE(again[0].features[i].type == rec.features[i].type);
        REQUIRE(again[0].features[i].parts == rec.features[i].parts);
        REQUIRE(again[0].features[i].strand == rec.features[i].strand);
        REQUIRE(again[0].features[i].qualifiers == rec.features[i].qualifiers);
    }
    // An empty description survives as empty: "." is the GenBank null value, not a
    // one-character description.
    NucRecord bare;
    bare.id = "bare";
    bare.sequence = "ACGTACGTACGT";
    const std::vector<NucRecord> bareAgain = readGenBank(writeGenBank(bare));
    REQUIRE(bareAgain.size() == 1);
    REQUIRE(bareAgain[0].description.empty());
    REQUIRE(bareAgain[0].sequence == bare.sequence);
}

TEST_CASE("FASTA reading and writing", "[bio][nuc]") {
    const std::string text = ">seq1 first record\nACGTACGT\nacgtn\n>seq2\nAUGCAU\n";
    const auto records = readNucFasta(text);
    REQUIRE(records.size() == 2);
    REQUIRE(records[0].id == "seq1");
    REQUIRE(records[0].description == "first record");
    REQUIRE(records[0].sequence == "ACGTACGTACGTN");
    REQUIRE(records[0].kind == NucKind::Dna);
    REQUIRE(records[1].kind == NucKind::Rna);
    REQUIRE(writeNucFasta(records[0]) == ">seq1 first record\nACGTACGTACGTN\n");
    // The sniffer routes both formats through one entry point.
    REQUIRE(readNucleic(text).size() == 2);
    REQUIRE(readNucleic(readAll(fixture("pUC19.gb")))[0].sequence.size() == 2686);
    REQUIRE_THROWS_AS(readNucleic("neither format"), Error);
}

TEST_CASE("six-frame translation and ORF finding, including a wrap", "[bio][nuc]") {
    const GeneticCodeTables tables = loadGeneticCodes(pack("genetic-codes.json"));
    const GeneticCode& code = *tables.find(1);
    const auto frames = sixFrameTranslation("ATGGCCATTGTAATGGGCCGCTGA", code);
    REQUIRE(frames.size() == 6);
    REQUIRE(frames[0] == "MAIVMGR*");
    REQUIRE(frames[3] == "SAAHYNGH");

    std::string circular;
    for (int i = 0; i < 30; ++i) circular += "GCT";
    circular += "TAA";
    circular += "ATG";
    OrfOptions options;
    options.minAminoAcids = 25;
    options.circular = true;
    const auto wrapped = findOrfs(circular, code, options);
    REQUIRE_FALSE(wrapped.empty());
    REQUIRE(wrapped[0].begin == static_cast<int>(circular.size()) - 3);
    // The end coordinate deliberately exceeds the length: that is how a wrap is
    // expressed without duplicating the sequence.
    REQUIRE(wrapped[0].end > static_cast<int>(circular.size()));
    REQUIRE(wrapped[0].protein.size() == 31);
    REQUIRE(wrapped[0].stopped);

    options.circular = false;
    REQUIRE(findOrfs(circular, code, options).empty());
}

TEST_CASE("the restriction pack is internally consistent", "[bio][nuc]") {
    const EnzymeSet set = loadEnzymes(pack("restriction-enzymes.json"));
    REQUIRE(set.all().size() >= 60);
    int blunt = 0, fivePrime = 0, threePrime = 0, ambiguous = 0, interrupted = 0;
    for (const auto& e : set.all()) {
        REQUIRE(e.topCut + e.bottomCut == static_cast<int>(e.recognition.size()));
        REQUIRE(reverseComplement(e.recognition) == e.recognition);
        if (e.overhang() == 0) ++blunt;
        if (e.overhang() > 0) ++fivePrime;
        if (e.overhang() < 0) ++threePrime;
        bool hasAmbiguity = false, hasN = false;
        for (char c : e.recognition) {
            if (c == 'N') hasN = true;
            else if (c != 'A' && c != 'C' && c != 'G' && c != 'T') hasAmbiguity = true;
        }
        if (hasAmbiguity) ++ambiguous;
        if (hasN) ++interrupted;
    }
    REQUIRE(blunt >= 1);
    REQUIRE(fivePrime >= 1);
    REQUIRE(threePrime >= 1);
    REQUIRE(ambiguous >= 3);
    REQUIRE(interrupted >= 1);
    // The catalogue notation, spot-checked against the three overhang classes.
    REQUIRE(set.find("EcoRI")->topCut == 1);
    REQUIRE(set.find("EcoRI")->overhang() == 4);
    REQUIRE(set.find("SmaI")->overhang() == 0);
    REQUIRE(set.find("PstI")->overhang() == -4);
    REQUIRE(set.find("EcoNI")->recognition == "CCTNNNNNAGG");
}

TEST_CASE("circular cut arithmetic is exact", "[bio][nuc]") {
    const EnzymeSet set = loadEnzymes(pack("restriction-enzymes.json"));
    const NucRecord& rec = puc19();

    const RestrictionDigest single = digest(rec, {"EcoRI"}, set);
    REQUIRE(single.sites.size() == 1);
    REQUIRE(single.sites[0].position == 396);
    REQUIRE(single.sites[0].overhang == 4);
    int sum = 0;
    for (int l : single.fragmentLengths) sum += l;
    REQUIRE(sum == 2686);
    REQUIRE(single.fragmentLengths.size() == 1);

    const RestrictionDigest pairDigest = digest(rec, {"EcoRI", "HindIII"}, set);
    REQUIRE(pairDigest.sites.size() == 2);
    sum = 0;
    for (int l : pairDigest.fragmentLengths) sum += l;
    REQUIRE(sum == 2686);
    REQUIRE(pairDigest.fragmentLengths.size() == 2);

    const RestrictionDigest triple = digest(rec, {"EcoRI", "HindIII", "PstI"}, set);
    sum = 0;
    for (int l : triple.fragmentLengths) sum += l;
    REQUIRE(sum == 2686);

    // An unknown enzyme is reported, not silently dropped.
    const RestrictionDigest unknown = digest(rec, {"EcoRI", "NotAnEnzyme"}, set);
    REQUIRE(unknown.sites.size() == 1);
    REQUIRE_FALSE(unknown.warnings.empty());

    // A site spanning the origin: found exactly once, and only when circular.
    NucRecord wrap;
    wrap.id = "origin-spanning";
    wrap.circular = true;
    wrap.sequence = "AATTCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCG";
    const int n = static_cast<int>(wrap.sequence.size());
    REQUIRE(findPattern(wrap.sequence, "GAATTC", true).size() == 1);
    REQUIRE(findPattern(wrap.sequence, "GAATTC", true)[0] == n - 1);
    REQUIRE(findPattern(wrap.sequence, "GAATTC", false).empty());
    const RestrictionDigest wrapped = digest(wrap, {"EcoRI"}, set);
    REQUIRE(wrapped.sites.size() == 1);
    sum = 0;
    for (int l : wrapped.fragmentLengths) sum += l;
    REQUIRE(sum == n);

    // A linear template of the same sequence has no site at all.
    NucRecord linear = wrap;
    linear.circular = false;
    REQUIRE(digest(linear, {"EcoRI"}, set).sites.empty());
    REQUIRE(digest(linear, {"EcoRI"}, set).fragmentLengths.size() == 1);
}
