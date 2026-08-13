#include "bio/Restriction.h"

#include <algorithm>
#include <fstream>
#include <mutex>

#include "bio/NucSeq.h"
#include "core/Assets.h"
#include "core/Error.h"

namespace biocad::bio {

namespace {

bool isPalindromic(const std::string& site) {
    const std::string rc = reverseComplement(site);
    return rc == site;
}

}  // namespace

const RestrictionEnzyme* EnzymeSet::find(std::string_view name) const {
    for (const auto& e : enzymes_) {
        if (e.name == name) return &e;
    }
    return nullptr;
}

std::vector<std::string> EnzymeSet::names() const {
    std::vector<std::string> out;
    out.reserve(enzymes_.size());
    for (const auto& e : enzymes_) out.push_back(e.name);
    return out;
}

EnzymeSet parseEnzymes(const nlohmann::json& j) {
    const int version = j.value("schemaVersion", 0);
    if (version != 1) {
        throw Error::parse("restriction-enzymes pack: unsupported schemaVersion " +
                           std::to_string(version));
    }
    const auto list = j.find("enzymes");
    if (list == j.end() || !list->is_array() || list->empty()) {
        throw Error::parse("restriction-enzymes pack: no enzymes");
    }
    EnzymeSet out;
    for (const auto& e : *list) {
        RestrictionEnzyme enzyme;
        enzyme.name = e.value("name", std::string{});
        enzyme.recognition = e.value("recognition", std::string{});
        enzyme.topCut = e.value("topCut", -1);
        enzyme.bottomCut = e.value("bottomCut", -1);
        enzyme.notes = e.value("notes", std::string{});
        enzyme.source = e.value("source", std::string{});
        if (enzyme.name.empty() || enzyme.recognition.empty()) {
            throw Error::parse("restriction-enzymes pack: entry without a name or recognition site");
        }
        for (char c : enzyme.recognition) {
            if (!isIupac(c)) {
                throw Error::parse(enzyme.name + ": recognition sequence contains '" +
                                   std::string(1, c) + "', which is not an IUPAC symbol");
            }
        }
        if (!isPalindromic(enzyme.recognition)) {
            throw Error::parse(enzyme.name + ": recognition sequence " + enzyme.recognition +
                               " is not palindromic; this loader models palindromic sites only");
        }
        const int n = static_cast<int>(enzyme.recognition.size());
        if (enzyme.topCut < 0 || enzyme.bottomCut < 0 || enzyme.topCut + enzyme.bottomCut != n) {
            throw Error::parse(enzyme.name + ": topCut + bottomCut must equal the site length (" +
                               std::to_string(n) + "), got " + std::to_string(enzyme.topCut) +
                               " + " + std::to_string(enzyme.bottomCut));
        }
        if (out.find(enzyme.name)) {
            throw Error::parse("restriction-enzymes pack: duplicate enzyme " + enzyme.name);
        }
        out.add(std::move(enzyme));
    }
    return out;
}

EnzymeSet loadEnzymes(const std::filesystem::path& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw Error::io("cannot open restriction-enzymes pack: " + file.string());
    nlohmann::json j;
    try {
        in >> j;
    } catch (const nlohmann::json::exception& e) {
        throw Error::parse(file.string() + ": " + e.what());
    }
    try {
        return parseEnzymes(j);
    } catch (const Error& e) {
        throw Error::parse(file.string() + ": " + e.message);
    }
}

const EnzymeSet& builtinEnzymes() {
    static EnzymeSet set;
    static std::once_flag once;
    static std::string failure;
    std::call_once(once, [] {
        const auto dir = core::assetDir("packs/nucleic");
        if (dir.empty()) {
            failure = "asset root not found";
            return;
        }
        try {
            set = loadEnzymes(dir / "restriction-enzymes.json");
        } catch (const Error& e) {
            failure = e.message;
        }
    });
    if (set.empty()) throw Error::io("restriction enzymes unavailable: " + failure);
    return set;
}

std::vector<RestrictionSite> findSites(const NucRecord& record, const RestrictionEnzyme& enzyme) {
    std::vector<RestrictionSite> out;
    const int n = static_cast<int>(record.sequence.size());
    if (n == 0) return out;
    for (int start : findPattern(record.sequence, enzyme.recognition, record.circular)) {
        RestrictionSite site;
        site.enzyme = enzyme.name;
        site.recognition = enzyme.recognition;
        site.overhang = enzyme.overhang();
        // The cut index is the fragment boundary: the phosphodiester bond broken
        // lies between base position-1 and base position. On a circular template
        // the boundary wraps with the site, which is why this is a modulo and not
        // a clamp.
        int cut = start + enzyme.topCut;
        if (record.circular) {
            cut %= n;
        } else if (cut <= 0 || cut >= n) {
            // A linear template cut exactly at an end does not separate anything.
            continue;
        }
        site.position = cut;
        site.strand = Strand::Forward;   // every enzyme in the pack is palindromic
        out.push_back(site);
    }
    return out;
}

std::vector<std::string> enzymesCutting(const NucRecord& record, int count, const EnzymeSet& set) {
    std::vector<std::string> out;
    for (const auto& e : set.all()) {
        if (static_cast<int>(findSites(record, e).size()) == count) out.push_back(e.name);
    }
    return out;
}

RestrictionDigest digest(const NucRecord& record, const std::vector<std::string>& enzymeNames) {
    return digest(record, enzymeNames, builtinEnzymes());
}

RestrictionDigest digest(const NucRecord& record, const std::vector<std::string>& enzymeNames,
                         const EnzymeSet& set) {
    RestrictionDigest out;
    out.recordId = record.id;
    out.circular = record.circular;
    const int n = static_cast<int>(record.sequence.size());

    bool ambiguous = false;
    for (char c : record.sequence) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            ambiguous = true;
            break;
        }
    }
    if (ambiguous) {
        out.warnings.push_back("the template contains IUPAC ambiguity codes; a site is reported "
                               "only where every base is consistent with the recognition "
                               "sequence, so real sites may be missed");
    }

    for (const auto& name : enzymeNames) {
        const RestrictionEnzyme* enzyme = set.find(name);
        if (!enzyme) {
            out.warnings.push_back("unknown enzyme '" + name + "': not in the pack, so it "
                                   "contributed no cuts to this digest");
            continue;
        }
        for (auto& site : findSites(record, *enzyme)) out.sites.push_back(std::move(site));
    }
    std::sort(out.sites.begin(), out.sites.end(),
              [](const RestrictionSite& a, const RestrictionSite& b) {
                  if (a.position != b.position) return a.position < b.position;
                  return a.enzyme < b.enzyme;
              });

    // Fragment boundaries are the distinct top-strand cut positions. Two enzymes
    // cutting the same bond produce one boundary, not two zero-length fragments.
    std::vector<int> cuts;
    for (const auto& s : out.sites) {
        if (cuts.empty() || cuts.back() != s.position) cuts.push_back(s.position);
    }
    std::vector<int> lengths;
    if (n > 0) {
        if (record.circular) {
            if (cuts.empty()) {
                lengths.push_back(n);   // uncut circle: one fragment of the full length
            } else {
                for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
                    lengths.push_back(cuts[i + 1] - cuts[i]);
                }
                // The last fragment wraps the origin; this term is what makes the
                // sum exactly n instead of n minus the first cut position.
                lengths.push_back(n - cuts.back() + cuts.front());
            }
        } else {
            int previous = 0;
            for (int cut : cuts) {
                lengths.push_back(cut - previous);
                previous = cut;
            }
            lengths.push_back(n - previous);
        }
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<int>());
    out.fragmentLengths = std::move(lengths);

    int sum = 0;
    for (int l : out.fragmentLengths) sum += l;
    if (sum != n) {
        // Unreachable by construction; kept because a silent off-by-one here is
        // exactly the bug this module exists to avoid.
        throw Error::internal("digest: fragment lengths sum to " + std::to_string(sum) +
                              " but the template is " + std::to_string(n) + " bp");
    }
    if (record.circular && out.sites.size() == 1) {
        out.warnings.push_back("a single cut in a circular template linearises it: the one "
                               "fragment is the full length");
    }
    return out;
}

}  // namespace biocad::bio
