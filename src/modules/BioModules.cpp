#include "modules/BioModules.h"

namespace biocad {

SequenceAlignment RealSequence::alignGlobal(const std::string& a, const std::string& b) const {
    return bioadapt::alignGlobal(a, b);
}

SequenceAlignment RealSequence::alignLocal(const std::string& a, const std::string& b) const {
    return bioadapt::alignLocal(a, b);
}

std::optional<bio::Structure> RealStructure::load(const std::filesystem::path& file) const {
    return bioadapt::loadStructureFile(file);
}

StructureComparison RealStructure::compare(const bio::Structure& ref,
                                           const bio::Structure& model) const {
    return bioadapt::compare(ref, model);
}

Quantity RealStructure::sasa(const bio::Structure& s) const {
    return bioadapt::sasa(s);
}

}  // namespace biocad
