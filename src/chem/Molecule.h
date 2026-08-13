// chem/Molecule.h - lightweight molecular graph (atoms + bonds) and element data.
// This is the in-house cheminformatics core (no RDKit dependency): a real SMILES
// parser feeds this graph, and Descriptors.* compute properties from it.
#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace biocad::chem {

struct ElementInfo {
    int z;
    const char* symbol;
    double mass;  // standard atomic weight
};

const ElementInfo* findElementBySymbol(std::string_view symbol);
const ElementInfo* elementByZ(int z);
double atomicMass(int z);
const char* symbolByZ(int z);

struct Atom {
    int  z = 6;             // atomic number
    int  charge = 0;
    int  explicitH = 0;     // hydrogens declared inside [ ]
    bool bracket = false;   // came from a bracket atom [ ... ]
    bool aromatic = false;
    int  implicitH = 0;     // computed for organic-subset atoms
    bool inRing = false;
    std::vector<int> nbr;    // neighbor atom indices (heavy atoms only)
    std::vector<int> bonds;  // incident bond indices

    [[nodiscard]] int totalH() const { return bracket ? explicitH : implicitH; }
    [[nodiscard]] int degree() const { return static_cast<int>(nbr.size()); }
};

struct Bond {
    int  a = 0;
    int  b = 0;
    double order = 1.0;     // 1, 2, 3; aromatic = 1.5
    bool aromatic = false;
    bool inRing = false;
    [[nodiscard]] int other(int x) const { return x == a ? b : a; }
};

struct Molecule {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    [[nodiscard]] int heavyAtomCount() const { return static_cast<int>(atoms.size()); }
    [[nodiscard]] bool empty() const { return atoms.empty(); }
};

}  // namespace biocad::chem
