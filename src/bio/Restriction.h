#pragma once

// Restriction mapping against a pack of type II endonucleases.
//
// The cut arithmetic is the whole content of this file, and it is where these
// tools get it wrong. Two invariants are asserted rather than assumed:
//
//   * the fragment lengths sum to the template length exactly, for a circular
//     template as well as a linear one; and
//   * a recognition site that spans the origin of a circular template is found
//     exactly once - not twice, which is what naive sequence doubling gives, and
//     not zero, which is what a linear search gives.
//
// Enzymes that cut asymmetrically outside their recognition sequence (BsaI, BsmBI
// and the rest of the type IIS set) are deliberately absent from the pack rather
// than approximated with a palindromic model.

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/Nucleic.h"

namespace biocad::bio {

struct RestrictionEnzyme {
    std::string name;
    std::string recognition;   // IUPAC, palindromic
    int         topCut = 0;    // bases between the site's 5' end and the top-strand cut
    int         bottomCut = 0; // same distance, in top-strand coordinates, for the bottom cut
    std::string notes;
    std::string source;

    // >0 is a 5' overhang, <0 a 3' overhang, 0 blunt.
    int overhang() const { return bottomCut - topCut; }
};

class EnzymeSet {
public:
    const RestrictionEnzyme* find(std::string_view name) const;
    std::vector<std::string> names() const;
    void add(RestrictionEnzyme e) { enzymes_.push_back(std::move(e)); }
    const std::vector<RestrictionEnzyme>& all() const { return enzymes_; }
    bool empty() const { return enzymes_.empty(); }

private:
    std::vector<RestrictionEnzyme> enzymes_;
};

// Rejects a non-palindromic recognition sequence or a cut that is not symmetric
// about the site centre, because both would mean the single-pattern search below
// is the wrong model for that enzyme.
EnzymeSet parseEnzymes(const nlohmann::json& j);
EnzymeSet loadEnzymes(const std::filesystem::path& file);
const EnzymeSet& builtinEnzymes();

// Digest with the named enzymes. An unknown name is a warning on the report, not
// a silent omission. `fragmentLengths` is descending, as a gel reads.
RestrictionDigest digest(const NucRecord& record, const std::vector<std::string>& enzymeNames);
RestrictionDigest digest(const NucRecord& record, const std::vector<std::string>& enzymeNames,
                         const EnzymeSet& set);

// Every site of one enzyme, without the fragment arithmetic.
std::vector<RestrictionSite> findSites(const NucRecord& record, const RestrictionEnzyme& enzyme);

// Enzymes from `set` that cut the record exactly `count` times - the query behind
// "what is unique in this vector".
std::vector<std::string> enzymesCutting(const NucRecord& record, int count, const EnzymeSet& set);

}  // namespace biocad::bio
