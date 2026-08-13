#include "data/Domain.h"

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

}  // namespace biocad
