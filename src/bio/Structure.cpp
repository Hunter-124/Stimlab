#include "bio/Structure.h"

#include <algorithm>
#include <cctype>
#include <string_view>
#include <unordered_map>

namespace biocad::bio {
namespace {

// The polymer residue table. MSE (selenomethionine) is deliberately here: it is a HETATM in
// every PDB entry that contains it, yet it is a genuine chain residue, and dropping it puts a
// hole in the sequence that then mis-registers every alignment downstream.
const std::unordered_map<std::string_view, char>& polymerTable() {
    static const std::unordered_map<std::string_view, char> kTable = {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
        {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
        {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
        {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
        {"MSE", 'M'}, {"SEC", 'U'}, {"PYL", 'O'},
        // Ambiguity codes that appear in SEQRES and in some depositions.
        {"ASX", 'B'}, {"GLX", 'Z'}, {"XLE", 'J'}, {"UNK", 'X'},
        // Nucleic acid polymer residues, so a DNA/RNA chain is not silently empty.
        {"DA", 'A'},  {"DC", 'C'},  {"DG", 'G'},  {"DT", 'T'},  {"DU", 'U'},
        {"A", 'A'},   {"C", 'C'},   {"G", 'G'},   {"U", 'U'},   {"I", 'I'},
    };
    return kTable;
}

std::string trimmed(const std::string& s) {
    const auto b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t");
    return s.substr(b, e - b + 1);
}

}  // namespace

const Atom* Residue::atom(const std::string& nameIn) const {
    // Atom names are stored verbatim (" CA ") because the fixed-column padding encodes the
    // element position, but callers naturally ask for "CA". Compare on the trimmed form.
    const std::string want = trimmed(nameIn);
    for (const Atom& a : atoms) {
        if (trimmed(a.name) == want) return &a;
    }
    return nullptr;
}

char Residue::oneLetter() const {
    std::string key = trimmed(name);
    for (char& c : key) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    const auto& table = polymerTable();
    const auto it = table.find(std::string_view{key});
    return it == table.end() ? 'X' : it->second;
}

const Model* Structure::model(int number) const {
    for (const Model& m : models) {
        if (m.number == number) return &m;
    }
    return nullptr;
}

std::size_t Structure::atomCount() const {
    std::size_t n = 0;
    for (const Model& m : models) {
        for (const Chain& c : m.chains) {
            for (const Residue& r : c.residues) n += r.atoms.size();
        }
    }
    return n;
}

std::vector<char> sequenceOf(const Chain& c) {
    std::vector<char> out;
    out.reserve(c.residues.size());
    for (const Residue& r : c.residues) {
        // Solvent and ligands are not sequence. A residue absent from the polymer table is
        // skipped entirely rather than emitted as 'X', because an 'X' would be interpreted
        // downstream as an unknown *polymer* position and shift the numbering.
        std::string key = trimmed(r.name);
        for (char& ch : key) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        const auto& table = polymerTable();
        const auto it = table.find(std::string_view{key});
        if (it == table.end()) continue;
        out.push_back(it->second);
    }
    return out;
}

}  // namespace biocad::bio
