#include "bio/NucIo.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <sstream>

#include "bio/NucSeq.h"
#include "core/Error.h"

namespace biocad::bio {

namespace {

std::string trim(std::string_view s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return std::string(s.substr(a, b - a));
}

std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> lines;
    std::string cur;
    for (char c : text) {
        if (c == '\n') {
            if (!cur.empty() && cur.back() == '\r') cur.pop_back();
            lines.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        if (cur.back() == '\r') cur.pop_back();
        lines.push_back(cur);
    }
    return lines;
}

// Splits a comma-separated argument list at depth 0 only, so join(complement(1..3),5..9)
// yields two arguments and not four.
std::vector<std::string> splitTopLevel(const std::string& s) {
    std::vector<std::string> out;
    int depth = 0;
    std::string cur;
    for (char c : s) {
        if (c == '(') ++depth;
        if (c == ')') --depth;
        if (c == ',' && depth == 0) {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!trim(cur).empty()) out.push_back(trim(cur));
    return out;
}

bool startsWithCall(const std::string& s, const char* name, std::string& inner) {
    const std::string prefix = std::string(name) + "(";
    if (s.size() <= prefix.size() || s.compare(0, prefix.size(), prefix) != 0) return false;
    if (s.back() != ')') return false;
    inner = s.substr(prefix.size(), s.size() - prefix.size() - 1);
    return true;
}

void parseLocationInto(const std::string& raw, ParsedLocation& out, bool complemented,
                       std::vector<std::string>* warnings) {
    const std::string s = trim(raw);
    if (s.empty()) return;
    std::string inner;
    if (startsWithCall(s, "complement", inner)) {
        parseLocationInto(inner, out, !complemented, warnings);
        return;
    }
    if (startsWithCall(s, "join", inner) || startsWithCall(s, "order", inner)) {
        if (s.compare(0, 5, "order") == 0) out.ordered = true;
        for (const auto& part : splitTopLevel(inner)) {
            parseLocationInto(part, out, complemented, warnings);
        }
        return;
    }
    if (complemented) out.strand = Strand::Reverse;

    std::string span = s;
    // A remote location names another entry; the bases are not in this record, so
    // it is recorded and skipped rather than resolved against the wrong sequence.
    const auto colon = span.find(':');
    if (colon != std::string::npos) {
        out.remoteAccession = span.substr(0, colon);
        if (warnings) {
            warnings->push_back("location references remote entry '" + out.remoteAccession +
                                "'; those bases are not in this record and were skipped");
        }
        return;
    }
    // A site between two bases, e.g. 1023^1024. Stored as a zero-length interval
    // so downstream arithmetic cannot mistake it for a base.
    const auto caret = span.find('^');
    if (caret != std::string::npos) {
        const int a = std::atoi(span.substr(0, caret).c_str());
        out.betweenBases = true;
        out.parts.emplace_back(a, a);
        if (warnings) warnings->push_back("location " + span + " is a site between bases");
        return;
    }
    int begin = 0, end = 0;
    const auto dots = span.find("..");
    std::string first = dots == std::string::npos ? span : span.substr(0, dots);
    std::string second = dots == std::string::npos ? span : span.substr(dots + 2);
    if (!first.empty() && first.front() == '<') {
        out.partial5 = true;
        first.erase(first.begin());
    }
    if (!first.empty() && first.front() == '>') {
        out.partial3 = true;
        first.erase(first.begin());
    }
    if (!second.empty() && second.back() == '>') {
        out.partial3 = true;
        second.pop_back();
    }
    if (!second.empty() && second.front() == '<') {
        out.partial5 = true;
        second.erase(second.begin());
    }
    if (!second.empty() && second.front() == '>') {
        out.partial3 = true;
        second.erase(second.begin());
    }
    begin = std::atoi(first.c_str());
    end = std::atoi(second.c_str());
    if (begin <= 0 || end <= 0) {
        throw Error::parse("unparsable feature location: " + s);
    }
    out.parts.emplace_back(begin - 1, end);   // GenBank is 1-based inclusive
}

}  // namespace

ParsedLocation parseLocation(const std::string& text, std::vector<std::string>* warnings) {
    ParsedLocation out;
    parseLocationInto(text, out, false, warnings);
    if (out.strand == Strand::Reverse) {
        // Both spellings of a reverse-strand spliced feature -
        // complement(join(a,b)) and join(complement(b),complement(a)) - denote the
        // same sequence. Storing parts ascending and complementing once, after
        // concatenation, collapses the two spellings onto one representation.
        std::sort(out.parts.begin(), out.parts.end());
    }
    if ((out.partial5 || out.partial3) && warnings) {
        warnings->push_back("location " + trim(text) +
                            " is partial; the stored bounds are the observed ends, not the "
                            "biological ones");
    }
    return out;
}

std::string formatLocation(const std::vector<std::pair<int, int>>& parts, Strand strand,
                           bool ordered) {
    if (parts.empty()) return {};
    std::string body;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i) body += ",";
        const int begin = parts[i].first + 1;
        const int end = parts[i].second;
        body += begin == end ? std::to_string(begin)
                             : std::to_string(begin) + ".." + std::to_string(end);
    }
    if (parts.size() > 1) body = (ordered ? "order(" : "join(") + body + ")";
    if (strand == Strand::Reverse) body = "complement(" + body + ")";
    return body;
}

// ------------------------------------------------------------------- FASTA

std::vector<NucRecord> readNucFasta(const std::string& text) {
    std::vector<NucRecord> out;
    NucRecord cur;
    bool open = false;
    std::string raw;
    const auto flush = [&] {
        if (!open) return;
        cur.sequence = sanitizeSequence(raw, &cur.warnings);
        // U in the payload is the only evidence a FASTA file gives about chemistry.
        cur.kind = raw.find('U') != std::string::npos || raw.find('u') != std::string::npos
                       ? NucKind::Rna
                       : NucKind::Dna;
        out.push_back(cur);
        cur = NucRecord{};
        raw.clear();
        open = false;
    };
    for (const auto& line : splitLines(text)) {
        if (!line.empty() && line[0] == '>') {
            flush();
            const std::string header = trim(line.substr(1));
            const auto sp = header.find_first_of(" \t");
            cur.id = sp == std::string::npos ? header : header.substr(0, sp);
            cur.description = sp == std::string::npos ? std::string{} : trim(header.substr(sp));
            open = true;
        } else if (open) {
            raw += line;
        }
    }
    flush();
    return out;
}

std::string writeNucFasta(const NucRecord& record, std::size_t lineWidth) {
    if (lineWidth == 0) throw Error::invalidArgument("writeNucFasta: lineWidth must be > 0");
    std::string out = ">" + record.id;
    if (!record.description.empty()) out += " " + record.description;
    out += "\n";
    for (std::size_t i = 0; i < record.sequence.size(); i += lineWidth) {
        out += record.sequence.substr(i, lineWidth);
        out += "\n";
    }
    return out;
}

// ----------------------------------------------------------------- GenBank

namespace {

// GenBank column layout for the feature table: the key occupies columns 6-20 and
// both locations and qualifiers begin at column 22. A line that starts at column
// 22 without a '/' continues the previous location or qualifier value.
constexpr std::size_t kFeatureKeyColumn = 5;
constexpr std::size_t kFeatureValueColumn = 21;

bool isFeatureStart(const std::string& line) {
    if (line.size() <= kFeatureKeyColumn) return false;
    for (std::size_t i = 0; i < kFeatureKeyColumn; ++i) {
        if (line[i] != ' ') return false;
    }
    return !std::isspace(static_cast<unsigned char>(line[kFeatureKeyColumn]));
}

bool isContinuation(const std::string& line) {
    if (line.size() <= kFeatureValueColumn) return false;
    for (std::size_t i = 0; i < kFeatureValueColumn; ++i) {
        if (line[i] != ' ') return false;
    }
    return true;
}

void appendQualifier(NucFeature& feature, const std::string& raw) {
    std::string body = raw;
    if (!body.empty() && body.front() == '/') body.erase(body.begin());
    const auto eq = body.find('=');
    std::string key = eq == std::string::npos ? body : body.substr(0, eq);
    std::string value = eq == std::string::npos ? std::string{} : body.substr(eq + 1);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        value = value.substr(1, value.size() - 2);
    }
    feature.qualifiers.emplace_back(trim(key), value);
}

// Continuation lines of a wrapped qualifier are joined with a space for free
// text, but /translation is a single token that NCBI wraps mid-word, so joining
// it with a space would corrupt the protein. That is the one special case.
void joinQualifierLine(NucFeature& feature, const std::string& piece) {
    if (feature.qualifiers.empty()) return;
    auto& q = feature.qualifiers.back();
    if (q.first == "translation") {
        q.second += piece;
    } else {
        if (!q.second.empty()) q.second += " ";
        q.second += piece;
    }
    if (q.second.size() >= 2 && q.second.front() == '"' && q.second.back() == '"') {
        q.second = q.second.substr(1, q.second.size() - 2);
    }
}

bool qualifierIsOpen(const NucFeature& feature) {
    if (feature.qualifiers.empty()) return false;
    const auto& v = feature.qualifiers.back().second;
    // A value that still carries its opening quote without a closing one is mid-flight.
    return !v.empty() && v.front() == '"';
}

}  // namespace

std::vector<NucRecord> readGenBank(const std::string& text) {
    std::vector<NucRecord> out;
    const auto lines = splitLines(text);

    NucRecord rec;
    bool inFeatures = false;
    bool inOrigin = false;
    bool haveRecord = false;
    std::string rawSequence;
    std::string pendingLocation;
    int declaredLength = 0;
    std::string lastHeaderKeyword;

    const auto finishLocation = [&] {
        if (pendingLocation.empty() || rec.features.empty()) {
            pendingLocation.clear();
            return;
        }
        auto& feature = rec.features.back();
        std::vector<std::string> locWarnings;
        const ParsedLocation loc = parseLocation(pendingLocation, &locWarnings);
        feature.parts = loc.parts;
        feature.strand = loc.strand;
        // Partiality is a property of the biology, so it is recorded on the
        // feature rather than only warned about and lost.
        if (loc.partial5) feature.qualifiers.emplace_back("biocad_partial_5prime", "true");
        if (loc.partial3) feature.qualifiers.emplace_back("biocad_partial_3prime", "true");
        if (loc.ordered) feature.qualifiers.emplace_back("biocad_location_order", "true");
        if (!loc.remoteAccession.empty()) {
            feature.qualifiers.emplace_back("biocad_remote_location", loc.remoteAccession);
        }
        for (auto& w : locWarnings) {
            rec.warnings.push_back(feature.type + ": " + w);
        }
        pendingLocation.clear();
    };

    const auto flush = [&] {
        if (!haveRecord) return;
        finishLocation();
        rec.sequence = sanitizeSequence(rawSequence, &rec.warnings);
        if (declaredLength > 0 && static_cast<int>(rec.sequence.size()) != declaredLength) {
            rec.warnings.push_back("LOCUS declares " + std::to_string(declaredLength) +
                                   " bp but ORIGIN carries " +
                                   std::to_string(rec.sequence.size()));
        }
        out.push_back(rec);
        rec = NucRecord{};
        rawSequence.clear();
        inFeatures = inOrigin = haveRecord = false;
        declaredLength = 0;
    };

    for (const auto& line : lines) {
        if (line.compare(0, 2, "//") == 0) {
            flush();
            continue;
        }
        if (line.compare(0, 5, "LOCUS") == 0) {
            flush();
            haveRecord = true;
            std::istringstream in(line.substr(5));
            std::string token;
            std::vector<std::string> tokens;
            while (in >> token) tokens.push_back(token);
            if (!tokens.empty()) rec.id = tokens[0];
            for (std::size_t i = 1; i < tokens.size(); ++i) {
                if (tokens[i] == "bp" && i >= 1) declaredLength = std::atoi(tokens[i - 1].c_str());
                if (tokens[i] == "circular") rec.circular = true;
                if (tokens[i] == "RNA" || tokens[i] == "mRNA" || tokens[i] == "rRNA" ||
                    tokens[i] == "tRNA" || tokens[i] == "ss-RNA" || tokens[i] == "ds-RNA") {
                    rec.kind = NucKind::Rna;
                }
            }
            lastHeaderKeyword = "LOCUS";
            continue;
        }
        if (!haveRecord) continue;

        if (inOrigin) {
            rawSequence += line;
            continue;
        }
        if (line.compare(0, 6, "ORIGIN") == 0) {
            finishLocation();
            inFeatures = false;
            inOrigin = true;
            continue;
        }
        if (line.compare(0, 8, "FEATURES") == 0) {
            inFeatures = true;
            lastHeaderKeyword.clear();
            continue;
        }

        if (inFeatures) {
            if (isFeatureStart(line)) {
                finishLocation();
                NucFeature feature;
                const std::string keyAndLoc = line.substr(kFeatureKeyColumn);
                feature.type = trim(keyAndLoc.substr(
                    0, std::min<std::size_t>(keyAndLoc.size(), kFeatureValueColumn - kFeatureKeyColumn)));
                rec.features.push_back(std::move(feature));
                pendingLocation = line.size() > kFeatureValueColumn
                                      ? trim(line.substr(kFeatureValueColumn))
                                      : std::string{};
                continue;
            }
            if (isContinuation(line)) {
                const std::string body = trim(line.substr(kFeatureValueColumn));
                if (!body.empty() && body[0] == '/' && !qualifierIsOpen(rec.features.back())) {
                    finishLocation();
                    appendQualifier(rec.features.back(), body);
                } else if (!pendingLocation.empty()) {
                    pendingLocation += body;   // a wrapped join() continues without a space
                } else if (!rec.features.empty()) {
                    joinQualifierLine(rec.features.back(), body);
                }
                continue;
            }
            // A non-indented line ends the feature table even without ORIGIN.
            finishLocation();
            inFeatures = false;
        }

        // Header section. Keywords start in column 1; continuation lines are
        // indented, and DEFINITION in particular is routinely two or three lines.
        if (!line.empty() && !std::isspace(static_cast<unsigned char>(line[0]))) {
            std::istringstream in(line);
            std::string keyword;
            in >> keyword;
            lastHeaderKeyword = keyword;
            const std::string value = line.size() > 12 ? trim(line.substr(12)) : std::string{};
            if (keyword == "DEFINITION") {
                // "." is the GenBank null value, and the writer below emits it for
                // an empty description. Reading it back as a literal description
                // would break the record -> file -> record identity.
                rec.description = value == "." ? std::string{} : value;
            } else if (keyword == "ACCESSION") {
                if (rec.id.empty()) rec.id = value;
            }
        } else if (lastHeaderKeyword == "DEFINITION") {
            const std::string cont = trim(line);
            if (!cont.empty()) rec.description += " " + cont;
        }
    }
    flush();
    return out;
}

std::string writeGenBank(const NucRecord& record) {
    std::string out;
    char buf[128];
    std::snprintf(buf, sizeof buf, "LOCUS       %-16s %9zu bp    %-6s  %-8s SYN\n",
                  record.id.c_str(), record.sequence.size(),
                  record.kind == NucKind::Rna ? "RNA" : "DNA",
                  record.circular ? "circular" : "linear");
    out += buf;

    // DEFINITION wraps at 79 columns with a 12-column continuation indent, which
    // is what the reader above expects to see.
    const auto wrapped = [&](const std::string& keyword, const std::string& value,
                             std::size_t indent) {
        std::string head = keyword;
        head.resize(std::max<std::size_t>(keyword.size(), indent), ' ');
        std::string line = head;
        std::size_t col = line.size();
        std::istringstream in(value);
        std::string word;
        bool first = true;
        while (in >> word) {
            if (!first && col + 1 + word.size() > 79) {
                out += line + "\n";
                line.assign(indent, ' ');
                col = indent;
                first = true;
            }
            if (!first) {
                line += " ";
                ++col;
            }
            line += word;
            col += word.size();
            first = false;
        }
        out += line + "\n";
    };
    wrapped("DEFINITION", record.description.empty() ? "." : record.description, 12);
    wrapped("ACCESSION", record.id, 12);
    out += "FEATURES             Location/Qualifiers\n";
    for (const auto& feature : record.features) {
        std::string key = feature.type;
        key.resize(std::max<std::size_t>(key.size(), 16), ' ');
        out += "     " + key + formatLocation(feature.parts, feature.strand,
                                              std::any_of(feature.qualifiers.begin(),
                                                          feature.qualifiers.end(),
                                                          [](const auto& q) {
                                                              return q.first ==
                                                                     "biocad_location_order";
                                                          })) + "\n";
        for (const auto& q : feature.qualifiers) {
            // The biocad_* qualifiers are this reader's record of the location
            // syntax; re-emitting them as free qualifiers *and* in the location
            // would duplicate them on the next read.
            if (q.first.compare(0, 7, "biocad_") == 0) continue;
            std::string text = "/" + q.first;
            if (!q.second.empty()) {
                const bool numeric = std::all_of(q.second.begin(), q.second.end(), [](char c) {
                    return std::isdigit(static_cast<unsigned char>(c)) != 0;
                });
                text += "=" + (numeric ? q.second : "\"" + q.second + "\"");
            }
            // /translation is one token that NCBI wraps mid-word, and the reader
            // rejoins it without a separator, so it is chopped at a fixed width.
            // Every other qualifier is rejoined with a single space, so it must be
            // broken only at a space or the round trip would insert one.
            const std::size_t width = 79 - kFeatureValueColumn;
            const bool chop = q.first == "translation";
            std::size_t pos = 0;
            while (pos < text.size()) {
                std::size_t take = std::min(width, text.size() - pos);
                if (!chop && pos + take < text.size()) {
                    const auto space = text.rfind(' ', pos + take);
                    if (space != std::string::npos && space > pos) take = space - pos;
                }
                out += std::string(kFeatureValueColumn, ' ') + text.substr(pos, take) + "\n";
                pos += take;
                while (!chop && pos < text.size() && text[pos] == ' ') ++pos;
            }
        }
    }
    out += "ORIGIN      \n";
    for (std::size_t i = 0; i < record.sequence.size(); i += 60) {
        std::snprintf(buf, sizeof buf, "%9zu", i + 1);
        out += buf;
        for (std::size_t j = 0; j < 60 && i + j < record.sequence.size(); j += 10) {
            out += " ";
            for (std::size_t k = 0; k < 10 && i + j + k < record.sequence.size(); ++k) {
                out += static_cast<char>(
                    std::tolower(static_cast<unsigned char>(record.sequence[i + j + k])));
            }
        }
        out += "\n";
    }
    out += "//\n";
    return out;
}

std::vector<NucRecord> readNucleic(const std::string& text) {
    for (const auto& line : splitLines(text)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '>') return readNucFasta(text);
        if (t.compare(0, 5, "LOCUS") == 0) return readGenBank(text);
        break;
    }
    throw Error::parse("input is neither FASTA (leading '>') nor GenBank (leading 'LOCUS')");
}

}  // namespace biocad::bio
