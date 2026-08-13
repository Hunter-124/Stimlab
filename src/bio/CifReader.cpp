#include "bio/CifReader.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace biocad::bio {
namespace {

// A CIF value plus whether it arrived quoted. The distinction matters: unquoted '.' means
// "inapplicable" and unquoted '?' means "unknown", but a QUOTED "?" is the literal character,
// which really does occur (chemical component names, free-text descriptions).
struct Value {
    std::string text;
    bool quoted = false;

    [[nodiscard]] bool missing() const {
        return !quoted && (text.empty() || text == "." || text == "?");
    }
};

struct Category {
    std::vector<std::string> tags;                // keyword part after "_category."
    std::vector<std::vector<Value>> rows;
};

struct Document {
    std::string dataBlock;
    std::unordered_map<std::string, Category> categories;   // key: "_atom_site"
};

bool isSpace(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

std::string lowered(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// The tokeniser. STAR is small but has four traps: '#' comments run to end of line, single and
// double quotes terminate only when the closing quote is followed by whitespace (so
// "5'-end" inside a bare token is not a quote), a ';' in column 1 opens a multi-line text field
// that ends at the next ';' in column 1, and tags/keywords are case-insensitive.
class Tokenizer {
public:
    explicit Tokenizer(const std::string& text) : text_(text) {}

    bool next(Value& out) {
        while (pos_ < text_.size()) {
            const char c = text_[pos_];
            if (isSpace(c)) {
                atLineStart_ = (c == '\n');
                ++pos_;
                continue;
            }
            if (c == '#') {
                while (pos_ < text_.size() && text_[pos_] != '\n') ++pos_;
                continue;
            }
            break;
        }
        if (pos_ >= text_.size()) return false;

        const char c = text_[pos_];
        if (c == ';' && atLineStart_) {
            ++pos_;
            std::string body;
            bool lineStart = false;
            while (pos_ < text_.size()) {
                if (lineStart && text_[pos_] == ';') {
                    ++pos_;
                    break;
                }
                lineStart = (text_[pos_] == '\n');
                body.push_back(text_[pos_]);
                ++pos_;
            }
            // Strip the leading newline that follows the opening ';' and the trailing one that
            // precedes the closing ';': they are delimiters, not content.
            if (!body.empty() && body.front() == '\n') body.erase(body.begin());
            if (!body.empty() && body.back() == '\n') body.pop_back();
            out = Value{body, true};
            atLineStart_ = false;
            return true;
        }
        if (c == '\'' || c == '"') {
            const char quote = c;
            ++pos_;
            std::string body;
            while (pos_ < text_.size()) {
                if (text_[pos_] == quote &&
                    (pos_ + 1 >= text_.size() || isSpace(text_[pos_ + 1]))) {
                    ++pos_;
                    break;
                }
                body.push_back(text_[pos_]);
                ++pos_;
            }
            out = Value{body, true};
            atLineStart_ = false;
            return true;
        }
        std::string body;
        while (pos_ < text_.size() && !isSpace(text_[pos_])) {
            body.push_back(text_[pos_]);
            ++pos_;
        }
        out = Value{body, false};
        atLineStart_ = false;
        return true;
    }

private:
    const std::string& text_;
    std::size_t pos_ = 0;
    bool atLineStart_ = true;
};

bool splitTag(const std::string& token, std::string& category, std::string& keyword) {
    if (token.empty() || token[0] != '_') return false;
    const auto dot = token.find('.');
    if (dot == std::string::npos) {
        category = lowered(token);
        keyword.clear();
        return true;
    }
    category = lowered(token.substr(0, dot));
    keyword = lowered(token.substr(dot + 1));
    return true;
}

Document parseDocument(const std::string& text) {
    Document doc;
    Tokenizer tk(text);
    std::vector<Value> pending;
    Value tok;
    bool have = tk.next(tok);
    while (have) {
        const std::string low = lowered(tok.text);
        if (!tok.quoted && low.rfind("data_", 0) == 0) {
            doc.dataBlock = tok.text.substr(5);
            have = tk.next(tok);
            continue;
        }
        if (!tok.quoted && (low == "loop_")) {
            // Collect the tag list, then read values until a token that cannot be a value.
            std::string category;
            std::vector<std::string> keywords;
            have = tk.next(tok);
            while (have && !tok.quoted && !tok.text.empty() && tok.text[0] == '_') {
                std::string cat, key;
                splitTag(tok.text, cat, key);
                if (category.empty()) category = cat;
                keywords.push_back(key);
                have = tk.next(tok);
            }
            std::vector<Value> flat;
            while (have) {
                const std::string vlow = lowered(tok.text);
                const bool terminator = !tok.quoted &&
                                        (tok.text.rfind("_", 0) == 0 || vlow == "loop_" ||
                                         vlow == "stop_" || vlow.rfind("data_", 0) == 0 ||
                                         vlow.rfind("save_", 0) == 0);
                if (terminator) break;
                flat.push_back(tok);
                have = tk.next(tok);
            }
            if (!keywords.empty()) {
                Category& c = doc.categories[category];
                if (c.tags.empty()) c.tags = keywords;
                const std::size_t width = keywords.size();
                for (std::size_t i = 0; i + width <= flat.size(); i += width) {
                    c.rows.emplace_back(flat.begin() + static_cast<std::ptrdiff_t>(i),
                                        flat.begin() + static_cast<std::ptrdiff_t>(i + width));
                }
            }
            continue;
        }
        if (!tok.quoted && !tok.text.empty() && tok.text[0] == '_') {
            std::string cat, key;
            splitTag(tok.text, cat, key);
            Value v;
            if (!tk.next(v)) break;
            // A non-loop item is a one-row category; merging it into the same structure means
            // the consumers below never care whether the file used loop_ or not, which is the
            // whole point since a single-row category is legally written either way.
            Category& c = doc.categories[cat];
            if (c.rows.empty()) c.rows.emplace_back();
            c.tags.push_back(key);
            c.rows.front().push_back(v);
            have = tk.next(tok);
            continue;
        }
        have = tk.next(tok);
    }
    return doc;
}

// Column lookup BY NAME. Returns npos when the tag is absent, which is normal: minimal files
// omit auth_* entirely and older files omit pdbx_PDB_model_num.
std::size_t columnOf(const Category& c, const char* keyword) {
    for (std::size_t i = 0; i < c.tags.size(); ++i) {
        if (c.tags[i] == keyword) return i;
    }
    return static_cast<std::size_t>(-1);
}

const Value* cell(const Category& c, const std::vector<Value>& row, const char* keyword) {
    const std::size_t i = columnOf(c, keyword);
    if (i == static_cast<std::size_t>(-1) || i >= row.size()) return nullptr;
    return &row[i];
}

std::string textOr(const Category& c, const std::vector<Value>& row, const char* keyword,
                   const std::string& fallback = {}) {
    const Value* v = cell(c, row, keyword);
    if (v == nullptr || v->missing()) return fallback;
    return v->text;
}

bool numberOf(const Category& c, const std::vector<Value>& row, const char* keyword, double& out) {
    const Value* v = cell(c, row, keyword);
    if (v == nullptr || v->missing()) return false;
    char* stop = nullptr;
    const double d = std::strtod(v->text.c_str(), &stop);
    if (stop != v->text.c_str() + v->text.size()) return false;
    out = d;
    return true;
}

bool intOf(const Category& c, const std::vector<Value>& row, const char* keyword, int& out) {
    const Value* v = cell(c, row, keyword);
    if (v == nullptr || v->missing()) return false;
    const char* begin = v->text.data();
    const char* end = begin + v->text.size();
    const auto res = std::from_chars(begin, end, out);   // handles a negative auth_seq_id
    return res.ec == std::errc{} && res.ptr == end;
}

char charOf(const Category& c, const std::vector<Value>& row, const char* keyword) {
    const Value* v = cell(c, row, keyword);
    if (v == nullptr || v->missing() || v->text.empty()) return ' ';
    return v->text[0];
}

struct Builder {
    Structure st;

    Model& model(int number) {
        for (Model& m : st.models) {
            if (m.number == number) return m;
        }
        st.models.push_back(Model{number, {}});
        return st.models.back();
    }

    void addAtom(int modelNumber, const std::string& chainId, const std::string& resName,
                 int authSeq, int labelSeq, char iCode, const Atom& atom) {
        Model& m = model(modelNumber);
        Chain* chain = nullptr;
        for (Chain& c : m.chains) {
            if (c.id == chainId) chain = &c;
        }
        if (chain == nullptr) {
            m.chains.push_back(Chain{chainId, {}});
            chain = &m.chains.back();
        }
        Residue* residue = chain->residues.empty() ? nullptr : &chain->residues.back();
        // Identity is (name, author number, insertion code) exactly as in the PDB reader, so the
        // same entry read from either format produces the same residue partition.
        const bool sameResidue = residue != nullptr && residue->name == resName &&
                                 residue->authSeqId == authSeq && residue->insertionCode == iCode;
        if (!sameResidue) {
            chain->residues.push_back(Residue{resName, authSeq, labelSeq, iCode, {}});
            residue = &chain->residues.back();
        }
        residue->atoms.push_back(atom);
    }
};

std::string paddedAtomName(const std::string& name, const std::string& element) {
    // mmCIF stores the bare atom name; the PDB stores it in a four-column field whose alignment
    // encodes the element. Re-create that padding so an atom read from either format compares
    // equal, using the same rule the PDB reader relies on: a one-letter element starts in
    // column 14, a two-letter element in column 13.
    if (name.size() >= 4) return name;
    std::string out;
    if (element.size() >= 2) {
        out = name;
    } else {
        out = " " + name;
    }
    while (out.size() < 4) out.push_back(' ');
    return out;
}

}  // namespace

Structure readCif(const std::string& text, std::string id, std::string source,
                  Annotations* annotations) {
    Builder b;
    const Document doc = parseDocument(text);
    b.st.id = id.empty() ? doc.dataBlock : std::move(id);
    b.st.source = std::move(source);

    const auto atomIt = doc.categories.find("_atom_site");
    if (atomIt == doc.categories.end()) {
        b.st.warnings.push_back("mmCIF has no _atom_site category");
    } else {
        const Category& c = atomIt->second;
        std::size_t rowNo = 0;
        for (const std::vector<Value>& row : c.rows) {
            ++rowNo;
            Atom a;
            if (!numberOf(c, row, "cartn_x", a.x) || !numberOf(c, row, "cartn_y", a.y) ||
                !numberOf(c, row, "cartn_z", a.z)) {
                b.st.warnings.push_back("_atom_site row " + std::to_string(rowNo) +
                                        ": unparsable coordinates");
                continue;
            }
            const std::string name = textOr(c, row, "auth_atom_id",
                                            textOr(c, row, "label_atom_id"));
            std::string resName = textOr(c, row, "auth_comp_id", textOr(c, row, "label_comp_id"));
            if (name.empty() || resName.empty()) {
                b.st.warnings.push_back("_atom_site row " + std::to_string(rowNo) +
                                        ": missing atom or residue name");
                continue;
            }
            a.element = textOr(c, row, "type_symbol");
            for (char& ch : a.element) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            a.name = paddedAtomName(name, a.element);
            a.altLoc = charOf(c, row, "label_alt_id");
            if (!numberOf(c, row, "occupancy", a.occupancy)) a.occupancy = 1.0;
            if (!numberOf(c, row, "b_iso_or_equiv", a.bFactor)) a.bFactor = 0.0;
            a.hetatm = (textOr(c, row, "group_pdb", "ATOM") == "HETATM");

            std::string chainId = textOr(c, row, "auth_asym_id", textOr(c, row, "label_asym_id"));
            if (chainId.empty()) chainId = " ";
            int authSeq = 0;
            if (!intOf(c, row, "auth_seq_id", authSeq) && !intOf(c, row, "label_seq_id", authSeq)) {
                b.st.warnings.push_back("_atom_site row " + std::to_string(rowNo) +
                                        ": unparsable residue sequence number");
                continue;
            }
            int labelSeq = 0;
            intOf(c, row, "label_seq_id", labelSeq);   // 0 stays "unknown", e.g. for waters
            const char iCode = charOf(c, row, "pdbx_pdb_ins_code");
            int modelNumber = 1;
            intOf(c, row, "pdbx_pdb_model_num", modelNumber);
            b.addAtom(modelNumber, chainId, resName, authSeq, labelSeq, iCode, a);
        }
    }

    if (annotations != nullptr) {
        const auto polyIt = doc.categories.find("_entity_poly");
        if (polyIt != doc.categories.end()) {
            const Category& c = polyIt->second;
            for (const std::vector<Value>& row : c.rows) {
                const std::string strands = textOr(c, row, "pdbx_strand_id");
                std::string seq = textOr(c, row, "pdbx_seq_one_letter_code_can",
                                         textOr(c, row, "pdbx_seq_one_letter_code"));
                // The one-letter code is a semicolon text block and is wrapped across lines.
                seq.erase(std::remove_if(seq.begin(), seq.end(),
                                         [](char ch) { return isSpace(ch); }),
                          seq.end());
                std::size_t begin = 0;
                while (begin <= strands.size() && !strands.empty()) {
                    const std::size_t comma = strands.find(',', begin);
                    const std::string one = strands.substr(
                        begin, comma == std::string::npos ? std::string::npos : comma - begin);
                    if (!one.empty()) annotations->seqres.push_back(SeqResChain{one, {}, seq});
                    if (comma == std::string::npos) break;
                    begin = comma + 1;
                }
            }
        }
        const auto confIt = doc.categories.find("_struct_conf");
        if (confIt != doc.categories.end()) {
            const Category& c = confIt->second;
            for (const std::vector<Value>& row : c.rows) {
                const std::string type = textOr(c, row, "conf_type_id");
                if (type.rfind("HELX", 0) != 0) continue;   // TURN_P and STRN are not helices
                HelixRecord h;
                h.id = textOr(c, row, "id");
                h.chainId = textOr(c, row, "beg_auth_asym_id",
                                   textOr(c, row, "beg_label_asym_id"));
                if (!intOf(c, row, "beg_auth_seq_id", h.startSeqId) ||
                    !intOf(c, row, "end_auth_seq_id", h.endSeqId)) {
                    b.st.warnings.push_back("_struct_conf: unparsable helix range");
                    continue;
                }
                h.startInsertionCode = charOf(c, row, "pdbx_beg_pdb_ins_code");
                h.endInsertionCode = charOf(c, row, "pdbx_end_pdb_ins_code");
                annotations->helices.push_back(h);
            }
        }
        const auto sheetIt = doc.categories.find("_struct_sheet_range");
        if (sheetIt != doc.categories.end()) {
            const Category& c = sheetIt->second;
            for (const std::vector<Value>& row : c.rows) {
                StrandRecord s;
                s.sheetId = textOr(c, row, "sheet_id");
                s.chainId = textOr(c, row, "beg_auth_asym_id",
                                   textOr(c, row, "beg_label_asym_id"));
                if (!intOf(c, row, "beg_auth_seq_id", s.startSeqId) ||
                    !intOf(c, row, "end_auth_seq_id", s.endSeqId)) {
                    b.st.warnings.push_back("_struct_sheet_range: unparsable strand range");
                    continue;
                }
                s.startInsertionCode = charOf(c, row, "pdbx_beg_pdb_ins_code");
                s.endInsertionCode = charOf(c, row, "pdbx_end_pdb_ins_code");
                annotations->strands.push_back(s);
            }
        }
        const auto connIt = doc.categories.find("_struct_conn");
        if (connIt != doc.categories.end()) {
            const Category& c = connIt->second;
            for (const std::vector<Value>& row : c.rows) {
                if (lowered(textOr(c, row, "conn_type_id")) != "disulf") continue;
                DisulfideRecord d;
                d.chainId1 = textOr(c, row, "ptnr1_auth_asym_id",
                                    textOr(c, row, "ptnr1_label_asym_id"));
                d.chainId2 = textOr(c, row, "ptnr2_auth_asym_id",
                                    textOr(c, row, "ptnr2_label_asym_id"));
                if (!intOf(c, row, "ptnr1_auth_seq_id", d.seqId1) ||
                    !intOf(c, row, "ptnr2_auth_seq_id", d.seqId2)) {
                    b.st.warnings.push_back("_struct_conn: unparsable disulfide");
                    continue;
                }
                d.insertionCode1 = charOf(c, row, "pdbx_ptnr1_pdb_ins_code");
                d.insertionCode2 = charOf(c, row, "pdbx_ptnr2_pdb_ins_code");
                annotations->disulfides.push_back(d);
            }
        }
    }

    if (b.st.models.empty()) b.st.models.push_back(Model{1, {}});
    return std::move(b.st);
}

Structure readCifFile(const std::filesystem::path& path, Annotations* annotations) {
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
    Structure st = readCif(buf.str(), {}, path.string(), annotations);
    if (st.id.empty()) st.id = path.stem().string();   // the data_ block wins; filename is a fallback
    return st;
}

}  // namespace biocad::bio
