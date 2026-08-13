#include "bio/Fasta.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <vector>

namespace biocad::bio {
namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::vector<std::string> splitSpaces(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string tok;
    while (in >> tok) out.push_back(tok);
    return out;
}

std::vector<std::string> splitOn(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

// A UniProt trailing field is exactly two uppercase letters followed by '='. Testing the shape
// rather than a fixed list would misfire on a protein name containing something like "AB=",
// so the caller checks membership in the known key set as well.
bool keyToken(const std::string& tok, std::string& key, std::string& firstValueWord) {
    if (tok.size() < 3 || tok[2] != '=') return false;
    if (!std::isupper(static_cast<unsigned char>(tok[0])) ||
        !std::isupper(static_cast<unsigned char>(tok[1]))) {
        return false;
    }
    key = tok.substr(0, 2);
    firstValueWord = tok.substr(3);
    return true;
}

bool knownUniProtKey(const std::string& key) {
    return key == "OS" || key == "OX" || key == "GN" || key == "PE" || key == "SV";
}

std::string joinWords(const std::vector<std::string>& words) {
    std::string out;
    for (const std::string& w : words) {
        if (!out.empty()) out.push_back(' ');
        out += w;
    }
    return out;
}

void parseUniProt(FastaRecord& r, const std::string& rest) {
    // Parsed RIGHT TO LEFT. Splitting left-to-right on spaces is wrong twice over: the protein
    // name contains spaces ("Sodium-dependent dopamine transporter"), and so do the OS= and GN=
    // values ("OS=Homo sapiens"). GN= is omitted entirely when the gene is unknown, so a
    // positional reading of the trailing fields is also wrong. Walking from the right, every
    // KEY= token closes the field whose words have accumulated to its right, and whatever words
    // are left over at the far left are the protein name.
    const std::vector<std::string> words = splitSpaces(rest);
    std::vector<std::string> buffer;   // words seen so far, in left-to-right order
    for (std::size_t i = words.size(); i-- > 0;) {
        std::string key, head;
        if (keyToken(words[i], key, head) && knownUniProtKey(key)) {
            std::vector<std::string> value;
            if (!head.empty()) value.push_back(head);
            value.insert(value.end(), buffer.begin(), buffer.end());
            const std::string joined = joinWords(value);
            if (key == "OS") r.organism = joined;
            else if (key == "OX") r.taxonId = joined;
            else if (key == "GN") r.geneName = joined;
            else if (key == "PE") r.evidence = joined;
            else if (key == "SV") r.version = joined;
            buffer.clear();
            continue;
        }
        buffer.insert(buffer.begin(), words[i]);
    }
    r.proteinName = joinWords(buffer);
}

void parseRcsb(FastaRecord& r, const std::vector<std::string>& fields) {
    // ">4M48_1|Chains A, B[auth C]|Protein name|Homo sapiens (9606)"
    r.entryId = trim(fields[0]);
    if (fields.size() > 1) {
        std::string chains = trim(fields[1]);
        const auto space = chains.find(' ');
        if (space != std::string::npos) chains = chains.substr(space + 1);   // drop "Chain(s)"
        for (std::string piece : splitOn(chains, ',')) {
            piece = trim(piece);
            if (piece.empty()) continue;
            // "B[auth C]" means the entity chain is B but the author called it C, and every
            // paper and every PDB file uses the author id, so that is the one kept.
            const auto open = piece.find("[auth ");
            if (open != std::string::npos) {
                const auto close = piece.find(']', open);
                if (close != std::string::npos) {
                    piece = trim(piece.substr(open + 6, close - open - 6));
                }
            }
            if (!piece.empty()) r.chainIds.push_back(piece);
        }
    }
    if (fields.size() > 2) r.proteinName = trim(fields[2]);
    if (fields.size() > 3) r.organism = trim(fields[3]);
}

void parseHeader(FastaRecord& r) {
    const std::string h = trim(r.header);
    const std::vector<std::string> fields = splitOn(h, '|');
    if (fields.size() >= 3 && (fields[0] == "sp" || fields[0] == "tr")) {
        r.kind = FastaHeaderKind::UniProt;
        r.database = fields[0];
        r.accession = trim(fields[1]);
        // The third field is "EntryName Protein name OS=...": the entry name is one token.
        std::string third = fields[2];
        // A '|' inside the description would have over-split; rejoin the tail defensively.
        for (std::size_t i = 3; i < fields.size(); ++i) third += "|" + fields[i];
        const std::string t = trim(third);
        const auto space = t.find(' ');
        r.entryName = space == std::string::npos ? t : t.substr(0, space);
        parseUniProt(r, space == std::string::npos ? std::string{} : trim(t.substr(space + 1)));
        return;
    }
    if (fields.size() >= 2 && !fields[0].empty() &&
        fields[1].rfind("Chain", 0) == 0) {
        r.kind = FastaHeaderKind::Rcsb;
        parseRcsb(r, fields);
        return;
    }
    r.kind = FastaHeaderKind::Plain;
    const auto space = h.find(' ');
    r.entryId = space == std::string::npos ? h : h.substr(0, space);
    if (space != std::string::npos) r.proteinName = trim(h.substr(space + 1));
}

}  // namespace

bool isSequenceCharacter(char c) {
    const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (u >= 'A' && u <= 'Z') return true;   // includes B Z J X U O
    return u == '*' || u == '-';
}

std::vector<FastaRecord> readFasta(const std::string& text, std::vector<std::string>* warnings) {
    std::vector<FastaRecord> out;
    std::istringstream in(text);
    std::string line;
    std::size_t lineNo = 0;
    bool sawHeader = false;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == ';') continue;   // the original FASTA comment character
        if (line[0] == '>') {
            out.push_back(FastaRecord{});
            out.back().header = trim(line.substr(1));
            parseHeader(out.back());
            sawHeader = true;
            continue;
        }
        if (!sawHeader) {
            if (warnings != nullptr) {
                warnings->push_back("line " + std::to_string(lineNo) +
                                    ": sequence data before any '>' header, ignored");
            }
            continue;
        }
        for (char c : line) {
            if (c == ' ' || c == '\t') continue;   // some tools pad sequence blocks
            if (!isSequenceCharacter(c)) {
                if (warnings != nullptr) {
                    warnings->push_back(std::string("line ") + std::to_string(lineNo) +
                                        ": dropped invalid sequence character '" + c + "'");
                }
                continue;
            }
            out.back().sequence.push_back(
                static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
        }
    }
    return out;
}

std::vector<FastaRecord> readFastaFile(const std::filesystem::path& path,
                                       std::vector<std::string>* warnings) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        if (warnings != nullptr) warnings->push_back("cannot open file: " + path.string());
        return {};
    }
    std::ostringstream buf;
    buf << f.rdbuf();
    return readFasta(buf.str(), warnings);
}

std::string writeFasta(const std::vector<FastaRecord>& records, std::size_t lineWidth) {
    if (lineWidth == 0) lineWidth = 60;
    std::string out;
    for (const FastaRecord& r : records) {
        out.push_back('>');
        out += r.header;
        out.push_back('\n');
        for (std::size_t i = 0; i < r.sequence.size(); i += lineWidth) {
            out += r.sequence.substr(i, lineWidth);
            out.push_back('\n');
        }
    }
    return out;
}

}  // namespace biocad::bio
