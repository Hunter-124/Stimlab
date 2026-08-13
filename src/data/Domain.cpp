#include "data/Domain.h"

#include <stdexcept>
#include <utility>

namespace biocad {

const char* verdictLabel(Verdict v) {
    switch (v) {
        case Verdict::Info:   return "Info";
        case Verdict::Good:   return "Good";
        case Verdict::Warn:   return "Warning";
        case Verdict::Danger: return "Danger";
    }
    return "Info";
}

const char* provenanceLabel(Provenance p) {
    switch (p) {
        case Provenance::Measured:    return "measured";
        case Provenance::Predicted:   return "predicted";
        case Provenance::Model:       return "model";
        case Provenance::Heuristic:   return "heuristic";
        case Provenance::NotComputed: return "not computed";
    }
    return "not computed";
}

Quantity makeQuantity(double value, std::string unit, double error, Provenance p,
                      std::string source) {
    if (p == Provenance::Heuristic && !unit.empty()) {
        throw std::invalid_argument(
            "heuristic quantity must not carry a physical unit (got \"" + unit +
            "\" from " + source + ")");
    }
    Quantity q;
    q.value = value;
    q.unit = std::move(unit);
    q.error = error;
    q.provenance = p;
    q.source = std::move(source);
    return q;
}

Quantity notComputed(std::string missingPrerequisite) {
    Quantity q;
    q.provenance = Provenance::NotComputed;
    q.source = std::move(missingPrerequisite);
    return q;
}

}  // namespace biocad
