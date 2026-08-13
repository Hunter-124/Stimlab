// Structure and sequence I/O: the PDB reader, the mmCIF reader, and FASTA.
//
// The load-bearing test is the cross-format one: the same synthetic entry is committed in both
// formats, and the two readers must agree atom for atom. That is what catches a column
// off-by-one in the fixed-column reader or a positional (rather than tag-named) column lookup
// in the mmCIF reader, neither of which any single-format test can see.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <string>

#include "bio/CifReader.h"
#include "bio/Fasta.h"
#include "bio/PdbReader.h"

using namespace biocad::bio;

namespace {

std::filesystem::path fixture(const char* name) {
    return std::filesystem::path(BIOCAD_TEST_FIXTURES) / name;
}

std::string sequenceString(const Chain& c) {
    const std::vector<char> s = sequenceOf(c);
    return std::string(s.begin(), s.end());
}

}  // namespace

TEST_CASE("PDB and mmCIF fixtures describe the same entry", "[bio][io]") {
    Annotations pdbAnn;
    Annotations cifAnn;
    const Structure p = readPdbFile(fixture("4m48_fragment.pdb"), &pdbAnn);
    const Structure c = readCifFile(fixture("4m48_fragment.cif"), &cifAnn);

    REQUIRE(p.warnings.empty());
    REQUIRE(c.warnings.empty());
    REQUIRE(p.atomCount() == 24);
    REQUIRE(p.atomCount() == c.atomCount());

    const Model* pm = p.model(1);
    const Model* cm = c.model(1);
    REQUIRE(pm != nullptr);
    REQUIRE(cm != nullptr);
    REQUIRE(pm->chains.size() == 2);
    REQUIRE(cm->chains.size() == pm->chains.size());

    for (std::size_t ci = 0; ci < pm->chains.size(); ++ci) {
        const Chain& pc = pm->chains[ci];
        const Chain& cc = cm->chains[ci];
        REQUIRE(pc.id == cc.id);
        REQUIRE(pc.residues.size() == cc.residues.size());
        for (std::size_t ri = 0; ri < pc.residues.size(); ++ri) {
            const Residue& pr = pc.residues[ri];
            const Residue& cr = cc.residues[ri];
            REQUIRE(pr.name == cr.name);
            REQUIRE(pr.authSeqId == cr.authSeqId);
            REQUIRE(pr.insertionCode == cr.insertionCode);
            REQUIRE(pr.atoms.size() == cr.atoms.size());
            for (std::size_t ai = 0; ai < pr.atoms.size(); ++ai) {
                const Atom& pa = pr.atoms[ai];
                const Atom& ca = cr.atoms[ai];
                REQUIRE(pa.name == ca.name);
                REQUIRE(pa.element == ca.element);
                REQUIRE(pa.altLoc == ca.altLoc);
                REQUIRE(pa.hetatm == ca.hetatm);
                REQUIRE(std::fabs(pa.x - ca.x) < 1e-6);
                REQUIRE(std::fabs(pa.y - ca.y) < 1e-6);
                REQUIRE(std::fabs(pa.z - ca.z) < 1e-6);
                REQUIRE(std::fabs(pa.occupancy - ca.occupancy) < 1e-6);
                REQUIRE(std::fabs(pa.bFactor - ca.bFactor) < 1e-6);
            }
        }
    }

    // label_seq_id exists only in mmCIF; the PDB has no such field and must report 0.
    REQUIRE(cm->chains[0].residues[2].labelSeqId == 3);
    REQUIRE(pm->chains[0].residues[2].labelSeqId == 0);

    REQUIRE(pdbAnn.seqres.size() == 2);
    REQUIRE(cifAnn.seqres.size() == 2);
    REQUIRE(pdbAnn.seqres[0].oneLetterCode == "GAS");
    REQUIRE(cifAnn.seqres[0].oneLetterCode == "GAS");
    REQUIRE(pdbAnn.helices.size() == 1);
    REQUIRE(cifAnn.helices.size() == 1);
    REQUIRE(pdbAnn.helices[0].chainId == "A");
    REQUIRE(pdbAnn.helices[0].endInsertionCode == 'A');
    REQUIRE(cifAnn.helices[0].endInsertionCode == 'A');
}

TEST_CASE("PDB atom-name alignment distinguishes carbon alpha from calcium", "[bio][io]") {
    const Structure p = readPdbFile(fixture("4m48_fragment.pdb"));
    const Model* m = p.model(1);
    REQUIRE(m != nullptr);
    const Chain& a = m->chains[0];
    const Chain& b = m->chains[1];

    // " CA " starts in column 14, so it is a one-letter element: carbon.
    const Atom* carbonAlpha = a.residues[0].atom("CA");
    REQUIRE(carbonAlpha != nullptr);
    REQUIRE(carbonAlpha->element == "C");
    REQUIRE_FALSE(carbonAlpha->hetatm);

    // "CA  " starts in column 13, so it is a two-letter element: calcium.
    const Residue* ion = nullptr;
    for (const Residue& r : b.residues) {
        if (r.name == "CA") ion = &r;
    }
    REQUIRE(ion != nullptr);
    REQUIRE(ion->atoms.size() == 1);
    REQUIRE(ion->atoms[0].name == "CA  ");
    REQUIRE(ion->atoms[0].element == "CA");
    REQUIRE(ion->atoms[0].hetatm);

    // Solvent and ions are not sequence.
    REQUIRE(sequenceString(a) == "GAS");
    REQUIRE(sequenceString(b) == "C");
}

TEST_CASE("Residue identity keeps negative numbering, insertion codes and altLocs", "[bio][io]") {
    const Structure p = readPdbFile(fixture("4m48_fragment.pdb"));
    const Chain& a = p.model(1)->chains[0];
    REQUIRE(a.residues.size() == 3);

    REQUIRE(a.residues[0].name == "GLY");
    REQUIRE(a.residues[0].authSeqId == -1);        // negative resSeq is legal and round-trips

    // 100 and 100A are DIFFERENT residues: the insertion code is part of identity.
    REQUIRE(a.residues[1].authSeqId == 100);
    REQUIRE(a.residues[1].insertionCode == ' ');
    REQUIRE(a.residues[2].authSeqId == 100);
    REQUIRE(a.residues[2].insertionCode == 'A');
    REQUIRE(a.residues[1].name == "ALA");
    REQUIRE(a.residues[2].name == "SER");

    // Both alternate locations survive; collapsing them would invent a residue with two
    // side-chain oxygens in impossible positions.
    const Residue& ser = a.residues[2];
    REQUIRE(ser.atoms.size() == 7);
    REQUIRE(ser.atoms[5].altLoc == 'A');
    REQUIRE(ser.atoms[6].altLoc == 'B');
    REQUIRE(std::fabs(ser.atoms[5].occupancy - 0.60) < 1e-9);
    REQUIRE(std::fabs(ser.atoms[6].occupancy - 0.40) < 1e-9);
}

TEST_CASE("Malformed and truncated input warns instead of throwing", "[bio][io]") {
    const std::string junk =
        "ATOM      1  N   GLY A  -1      12.345  -3.210   8.001  1.00 21.30\n"
        "ATOM      2  CA  GLY A  -1      not-a-number  -2.005   8.440\n"
        "ATOM      3  C   GLY A  ZZ      14.560  -2.330   8.810\n"
        "HETATM    4\n"
        "ATOM      5  O   GLY A  -1      15.010\n"
        "ATOM      6";
    Structure t;
    REQUIRE_NOTHROW(t = readPdb(junk, "junk"));
    REQUIRE(t.atomCount() == 1);          // the one good line still loads
    REQUIRE(t.warnings.size() == 5);

    // An mmCIF loop cut off mid-row keeps the whole rows and drops the partial one.
    const std::string truncated =
        "data_X\nloop_\n_atom_site.group_PDB\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n"
        "_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "ATOM N N ALA A 1 1.0 2.0 3.0\n"
        "ATOM C CA ALA A 1 4.0 5.0\n";
    Structure c;
    REQUIRE_NOTHROW(c = readCif(truncated));
    REQUIRE(c.atomCount() == 1);
    REQUIRE(c.id == "X");
}

TEST_CASE("UniProt FASTA headers parse right to left", "[bio][io]") {
    const std::string text =
        ">sp|Q01959|SC6A3_HUMAN Sodium-dependent dopamine transporter OS=Homo sapiens "
        "OX=9606 GN=SLC6A3 PE=1 SV=2\n"
        "MSKSKCSVGLMSSVVAPAKEPNAVGPKE*\n"
        ">sp|P0DTC2|X_TEST Uncharacterized protein with a long name OS=Escherichia coli "
        "(strain K12) OX=83333 PE=4 SV=1\n"
        "BZJXUO-ACDEFG\n"
        ">4M48_1|Chains A, B[auth C]|Sodium-dependent dopamine transporter|Homo sapiens (9606)\n"
        "GASC\n";
    std::vector<std::string> warnings;
    const std::vector<FastaRecord> recs = readFasta(text, &warnings);
    REQUIRE(recs.size() == 3);
    REQUIRE(warnings.empty());

    REQUIRE(recs[0].kind == FastaHeaderKind::UniProt);
    REQUIRE(recs[0].accession == "Q01959");
    REQUIRE(recs[0].entryName == "SC6A3_HUMAN");
    REQUIRE(recs[0].proteinName == "Sodium-dependent dopamine transporter");
    REQUIRE(recs[0].organism == "Homo sapiens");
    REQUIRE(recs[0].geneName == "SLC6A3");
    REQUIRE(recs[0].sequence.back() == '*');       // the stop codon is retained

    // The hard case: a multi-word protein name, a multi-word organism, and NO GN= at all.
    REQUIRE(recs[1].proteinName == "Uncharacterized protein with a long name");
    REQUIRE(recs[1].organism == "Escherichia coli (strain K12)");
    REQUIRE(recs[1].geneName.empty());
    REQUIRE(recs[1].taxonId == "83333");
    REQUIRE(recs[1].evidence == "4");
    REQUIRE(recs[1].version == "1");
    REQUIRE(recs[1].sequence == "BZJXUO-ACDEFG");  // B Z J X U O plus the gap character

    REQUIRE(recs[2].kind == FastaHeaderKind::Rcsb);
    REQUIRE(recs[2].entryId == "4M48_1");
    REQUIRE(recs[2].chainIds.size() == 2);
    REQUIRE(recs[2].chainIds[0] == "A");
    REQUIRE(recs[2].chainIds[1] == "C");           // "[auth C]" wins, that is what papers cite

    const std::string round = writeFasta({recs[1]}, 5);
    const std::vector<FastaRecord> again = readFasta(round);
    REQUIRE(again.size() == 1);
    REQUIRE(again[0].sequence == recs[1].sequence);
    REQUIRE(again[0].proteinName == recs[1].proteinName);
    REQUIRE(again[0].organism == recs[1].organism);
}
