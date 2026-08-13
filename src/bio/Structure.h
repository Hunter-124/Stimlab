#pragma once

// The protein object model. A protein is a hierarchy (model -> chain -> residue -> atom)
// with identity rules that small-molecule code does not have: residues are identified by
// an author number PLUS an insertion code, coordinates may carry alternate locations, and
// two independent numbering schemes (author and mmCIF label) coexist in the same file.
//
// Deliberate omission: there is no constructor from chem::Conformer and there never may be.
// chem::Conformer is a distance-geometry embedding of a small molecule; its coordinates are
// a plausible 3D layout, not a measured structure. Feeding one into TM-score, lDDT or GDT
// would produce a number that looks like a structure comparison and means nothing. The
// separation is enforced by the absence of any such constructor or overload, so a misuse is
// a compile error rather than a plausible-looking result.

#include <cstddef>
#include <string>
#include <vector>

namespace biocad::bio {

struct Atom {
    std::string name;        // PDB atom name, e.g. " CA "
    std::string element;     // "C", "N", "SE"
    char        altLoc = ' ';
    double      x = 0, y = 0, z = 0;
    double      occupancy = 1.0;
    double      bFactor = 0.0;
    bool        hetatm = false;
};

struct Residue {
    std::string name;        // "ALA", "HOH"
    int         authSeqId = 0;    // author numbering (what papers cite)
    int         labelSeqId = 0;   // mmCIF label numbering; 0 when unknown
    char        insertionCode = ' ';  // part of residue IDENTITY, not decoration
    std::vector<Atom> atoms;

    [[nodiscard]] const Atom* atom(const std::string& name) const;
    [[nodiscard]] char oneLetter() const;   // 'X' for anything unknown
};

struct Chain {
    std::string id;
    std::vector<Residue> residues;
};

struct Model {
    int number = 1;
    std::vector<Chain> chains;
};

struct Structure {
    std::string id;                       // "4M48"
    std::string source;                   // file path or URL it came from
    std::vector<Model> models;
    std::vector<std::string> warnings;    // recoverable parse problems, shown in the UI

    [[nodiscard]] const Model* model(int number = 1) const;
    [[nodiscard]] std::size_t atomCount() const;

    // Residue numbering is ambiguous: any UI that shows a number must state which.
    enum class Numbering { Author, Label };
};

std::vector<char> sequenceOf(const Chain& c);   // one-letter, polymer residues only

}  // namespace biocad::bio
