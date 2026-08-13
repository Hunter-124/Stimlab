#include "bio/PdbReader.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>

namespace biocad::bio {
namespace {

// Columns in the PDB spec are 1-based and inclusive. Every accessor goes through this so an
// off-by-one lives in exactly one place, and a short (truncated) line yields a short field
// instead of undefined behaviour.
std::string_view cols(std::string_view line, std::size_t first, std::size_t last) {
    if (first == 0 || first > last) return {};
    if (line.size() < first) return {};
    const std::size_t begin = first - 1;
    const std::size_t end = std::min(line.size(), last);
    return line.substr(begin, end - begin);
}

std::string_view trim(std::string_view s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '\r')) s.remove_prefix(1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.remove_suffix(1);
    return s;
}

bool parseInt(std::string_view s, int& out) {
    const std::string_view t = trim(s);
    if (t.empty()) return false;
    // from_chars handles the leading '-': resSeq is legally NEGATIVE (expression tags and
    // pre-mature-protein numbering routinely start below 1), so a reader that treats '-' as
    // junk drops the first residues of a great many entries.
    const auto* begin = t.data();
    const auto* end = t.data() + t.size();
    const auto res = std::from_chars(begin, end, out);
    return res.ec == std::errc{} && res.ptr == end;
}

bool parseDouble(std::string_view s, double& out) {
    const std::string_view t = trim(s);
    if (t.empty()) return false;
    // libstdc++ 11+ has floating-point from_chars, but strtod keeps this readable and is not
    // hot: a 100k-atom file is three strtod calls per line.
    const std::string tmp(t);
    char* stop = nullptr;
    const double v = std::strtod(tmp.c_str(), &stop);
    if (stop != tmp.c_str() + tmp.size()) return false;
    out = v;
    return true;
}

char upper(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }

// Element inference from the atom-name field (columns 13-16) when columns 77-78 are blank.
//
// The rule is positional, not textual: a one-letter element's name is right-shifted so that it
// begins in column 14, while a two-letter element begins in column 13. That is the whole
// difference between " CA " (carbon alpha of an amino acid) and "CA  " (a calcium ion).
// Getting it wrong silently turns a protein backbone into a pile of metal, and no downstream
// check catches it because both are chemically valid atoms.
std::string elementFromName(std::string_view rawName) {
    if (rawName.size() < 2) return std::string(trim(rawName));
    const char c13 = rawName[0];
    const char c14 = rawName[1];
    // Hydrogens are frequently written "1HB " / "2HG1": a digit in column 13 is a count prefix,
    // not an element, so the element is the column-14 character.
    if (c13 == ' ' || std::isdigit(static_cast<unsigned char>(c13))) {
        return c14 == ' ' ? std::string{} : std::string(1, upper(c14));
    }
    if (std::isalpha(static_cast<unsigned char>(c14))) {
        return std::string{upper(c13), upper(c14)};
    }
    return std::string(1, upper(c13));
}

struct Builder {
    Structure st;
    Model* current = nullptr;
    Chain* chain = nullptr;
    Residue* residue = nullptr;
    std::string chainKey;
    std::string residueKey;

    void beginModel(int number) {
        st.models.push_back(Model{number, {}});
        current = &st.models.back();
        chain = nullptr;
        residue = nullptr;
        chainKey.clear();
        residueKey.clear();
    }

    Model& model() {
        if (current == nullptr) beginModel(1);
        return *current;
    }

    void addAtom(const std::string& chainId, const std::string& resName, int authSeq,
                 char iCode, const Atom& atom) {
        Model& m = model();
        if (chain == nullptr || chainKey != chainId) {
            // A chain id that reappears after other chains is the same chain continuing (HETATM
            // groups are commonly emitted after every TER), so reuse an existing Chain.
            Chain* found = nullptr;
            for (Chain& c : m.chains) {
                if (c.id == chainId) found = &c;
            }
            if (found == nullptr) {
                m.chains.push_back(Chain{chainId, {}});
                found = &m.chains.back();
            }
            chain = found;
            chainKey = chainId;
            residue = nullptr;
            residueKey.clear();
        }
        // Residue identity is (name, author number, insertion code). 100 and 100A are two
        // different residues and must never be merged; the insertion code is identity, not
        // decoration.
        std::string key = resName;
        key += '|';
        key += std::to_string(authSeq);
        key += '|';
        key += iCode;
        if (residue == nullptr || residueKey != key) {
            chain->residues.push_back(Residue{resName, authSeq, 0, iCode, {}});
            residue = &chain->residues.back();
            residueKey = key;
        }
        residue->atoms.push_back(atom);
    }

    void endChain() {
        // TER closes the polymer; drop the cursors so nothing after it merges into the last
        // residue by accident.
        chain = nullptr;
        residue = nullptr;
        chainKey.clear();
        residueKey.clear();
    }
};

void warn(Structure& st, std::size_t lineNo, std::string_view what, std::string_view line) {
    std::ostringstream os;
    os << "line " << lineNo << ": " << what << ": \"" << std::string(trim(line)) << '"';
    st.warnings.push_back(os.str());
}

bool rangeEndpoints(std::string_view line, std::size_t chainCol, std::size_t seqFirst,
                    std::size_t seqLast, std::size_t iCodeCol, std::string& chainOut,
                    int& seqOut, char& iCodeOut) {
    chainOut = std::string(trim(cols(line, chainCol, chainCol)));
    if (!parseInt(cols(line, seqFirst, seqLast), seqOut)) return false;
    const std::string_view ic = cols(line, iCodeCol, iCodeCol);
    iCodeOut = ic.empty() ? ' ' : ic[0];
    return true;
}

}  // namespace

Structure readPdb(const std::string& text, std::string id, std::string source,
                  Annotations* annotations) {
    Builder b;
    b.st.id = std::move(id);
    b.st.source = std::move(source);

    std::istringstream in(text);
    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string_view record = trim(cols(line, 1, 6));
        if (record.empty()) continue;

        if (record == "HEADER") {
            const std::string_view code = trim(cols(line, 63, 66));
            if (b.st.id.empty() && !code.empty()) b.st.id = std::string(code);
            continue;
        }
        if (record == "MODEL") {
            int number = 0;
            if (!parseInt(cols(line, 11, 14), number)) {
                warn(b.st, lineNo, "MODEL record has no serial number", line);
                number = static_cast<int>(b.st.models.size()) + 1;
            }
            b.beginModel(number);
            continue;
        }
        if (record == "ENDMDL") {
            b.endChain();
            continue;
        }
        if (record == "TER") {
            b.endChain();
            continue;
        }
        if (record == "ATOM" || record == "HETATM") {
            const std::string_view rawName = cols(line, 13, 16);
            const std::string_view resName = trim(cols(line, 18, 20));
            if (trim(rawName).empty() || resName.empty()) {
                warn(b.st, lineNo, "atom record missing atom or residue name", line);
                continue;
            }
            int authSeq = 0;
            if (!parseInt(cols(line, 23, 26), authSeq)) {
                warn(b.st, lineNo, "unparsable residue sequence number", line);
                continue;
            }
            Atom a;
            if (!parseDouble(cols(line, 31, 38), a.x) || !parseDouble(cols(line, 39, 46), a.y) ||
                !parseDouble(cols(line, 47, 54), a.z)) {
                warn(b.st, lineNo, "unparsable coordinates", line);
                continue;
            }
            // Occupancy and B-factor are optional in practice; a missing field is not a reason
            // to discard a perfectly good coordinate.
            if (!parseDouble(cols(line, 55, 60), a.occupancy)) a.occupancy = 1.0;
            if (!parseDouble(cols(line, 61, 66), a.bFactor)) a.bFactor = 0.0;
            // Columns 67-76 are UNDEFINED in PDB v3.3 and are never parsed. Older software wrote
            // footprints, segment ids and site labels there; reading them yields garbage.
            const std::string_view elementField = trim(cols(line, 77, 78));
            a.name = std::string(rawName);
            a.element = elementField.empty() ? elementFromName(rawName)
                                             : std::string(elementField);
            for (char& c : a.element) c = upper(c);
            // altLoc is retained verbatim. Collapsing it here would silently merge two mutually
            // exclusive conformations into one impossible residue; selecting an altLoc is a
            // decision for the caller, not the parser.
            const std::string_view alt = cols(line, 17, 17);
            a.altLoc = alt.empty() ? ' ' : alt[0];
            a.hetatm = (record == "HETATM");

            std::string chainId(trim(cols(line, 22, 22)));
            if (chainId.empty()) chainId = " ";   // a blank chain id is legal and is its own chain
            const std::string_view icField = cols(line, 27, 27);
            const char iCode = icField.empty() ? ' ' : icField[0];
            b.addAtom(chainId, std::string(resName), authSeq, iCode, a);
            continue;
        }
        if (record == "SEQRES") {
            if (annotations == nullptr) continue;
            std::string chainId(trim(cols(line, 12, 12)));
            if (chainId.empty()) chainId = " ";
            SeqResChain* target = nullptr;
            for (SeqResChain& s : annotations->seqres) {
                if (s.chainId == chainId) target = &s;
            }
            if (target == nullptr) {
                annotations->seqres.push_back(SeqResChain{chainId, {}, {}});
                target = &annotations->seqres.back();
            }
            // Thirteen residue slots of four columns each, starting at column 20.
            for (std::size_t i = 0; i < 13; ++i) {
                const std::size_t first = 20 + i * 4;
                const std::string_view name = trim(cols(line, first, first + 2));
                if (name.empty()) continue;
                target->residueNames.emplace_back(name);
                // Mirror the mmCIF one-letter form so the two readers are directly comparable.
                Residue probe;
                probe.name = std::string(name);
                target->oneLetterCode.push_back(probe.oneLetter());
            }
            continue;
        }
        if (record == "HELIX") {
            if (annotations == nullptr) continue;
            HelixRecord h;
            h.id = std::string(trim(cols(line, 12, 14)));
            if (!rangeEndpoints(line, 20, 22, 25, 26, h.chainId, h.startSeqId,
                                h.startInsertionCode) ||
                !rangeEndpoints(line, 32, 34, 37, 38, h.chainId, h.endSeqId, h.endInsertionCode)) {
                warn(b.st, lineNo, "unparsable HELIX range", line);
                continue;
            }
            annotations->helices.push_back(h);
            continue;
        }
        if (record == "SHEET") {
            if (annotations == nullptr) continue;
            StrandRecord s;
            s.sheetId = std::string(trim(cols(line, 12, 14)));
            if (!rangeEndpoints(line, 22, 23, 26, 27, s.chainId, s.startSeqId,
                                s.startInsertionCode) ||
                !rangeEndpoints(line, 33, 34, 37, 38, s.chainId, s.endSeqId, s.endInsertionCode)) {
                warn(b.st, lineNo, "unparsable SHEET range", line);
                continue;
            }
            annotations->strands.push_back(s);
            continue;
        }
        if (record == "SSBOND") {
            if (annotations == nullptr) continue;
            DisulfideRecord d;
            if (!rangeEndpoints(line, 16, 18, 21, 22, d.chainId1, d.seqId1, d.insertionCode1) ||
                !rangeEndpoints(line, 30, 32, 35, 36, d.chainId2, d.seqId2, d.insertionCode2)) {
                warn(b.st, lineNo, "unparsable SSBOND", line);
                continue;
            }
            annotations->disulfides.push_back(d);
            continue;
        }
        // Any other record type is metadata this reader does not model; ignoring it is correct
        // and must not produce a warning, or a normal entry would emit hundreds.
    }

    if (b.st.models.empty()) b.st.models.push_back(Model{1, {}});
    return std::move(b.st);
}

Structure readPdbFile(const std::filesystem::path& path, Annotations* annotations) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        Structure st;
        st.source = path.string();
        st.warnings.push_back("cannot open file: " + path.string());
        st.models.push_back(Model{1, {}});
        return st;
    }
    std::ostringstream buf;
    buf << f.rdbuf();
    // The HEADER record is the authoritative id; the filename is only a fallback, because
    // downloaded files are routinely renamed.
    Structure st = readPdb(buf.str(), {}, path.string(), annotations);
    if (st.id.empty()) st.id = path.stem().string();
    return st;
}

}  // namespace biocad::bio
