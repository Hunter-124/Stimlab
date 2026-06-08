#include "chem/Molecule.h"

#include <array>

namespace stimlab::chem {

namespace {

// Standard atomic weights (IUPAC) for the elements StimLab parses.
constexpr std::array<ElementInfo, 25> kElements = {{
    {1, "H", 1.008},    {5, "B", 10.81},    {6, "C", 12.011},   {7, "N", 14.007},
    {8, "O", 15.999},   {9, "F", 18.998},   {11, "Na", 22.990}, {12, "Mg", 24.305},
    {13, "Al", 26.982}, {14, "Si", 28.085}, {15, "P", 30.974},  {16, "S", 32.06},
    {17, "Cl", 35.45},  {19, "K", 39.098},  {20, "Ca", 40.078}, {26, "Fe", 55.845},
    {33, "As", 74.922}, {34, "Se", 78.971}, {35, "Br", 79.904}, {53, "I", 126.904},
    {3, "Li", 6.94},    {30, "Zn", 65.38},  {78, "Pt", 195.084},{2, "He", 4.003},
    {10, "Ne", 20.180},
}};

}  // namespace

const ElementInfo* findElementBySymbol(std::string_view symbol) {
    for (const auto& e : kElements)
        if (symbol == e.symbol) return &e;
    return nullptr;
}

const ElementInfo* elementByZ(int z) {
    for (const auto& e : kElements)
        if (e.z == z) return &e;
    return nullptr;
}

double atomicMass(int z) {
    const ElementInfo* e = elementByZ(z);
    return e ? e->mass : 0.0;
}

const char* symbolByZ(int z) {
    const ElementInfo* e = elementByZ(z);
    return e ? e->symbol : "?";
}

}  // namespace stimlab::chem
