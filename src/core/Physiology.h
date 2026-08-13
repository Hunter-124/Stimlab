// core/Physiology.h - the single source of BioCAD's physiological constants.
//
// WHY THIS EXISTS. Hepatic blood flow appeared as the literal 90.0 inside
// chem/AdmetModel.h's well-stirred model while the drug-interaction pack carried
// 97 L/h. A constant that appears in two places has already diverged, and the
// version with no citation is a guess wearing units. assets/packs/physiology.json
// is now the only place these numbers exist, and this loader is the only way to
// read them.
//
// It lives in core, not in packs, for the same reason core/Assets.h does:
// biocad_chem must be able to reach it, and biocad_contracts already depends on
// biocad_chem, so a dependency on biocad_packs would close a cycle.
//
// A MISSING PACK IS NOT A DEFAULT. When the file cannot be read, `loaded` is false
// and every flow is zero, so a consumer divides by nothing and must report the
// failure. Substituting a built-in fallback here would silently restore exactly the
// uncited literal this file exists to delete.
#pragma once

#include <map>
#include <string>

namespace biocad::core {

struct Physiology {
    bool        loaded = false;
    double      hepaticBloodFlowLPerH = 0.0;      // Qh
    double      enterocyteBloodFlowLPerH = 0.0;   // Qen
    double      glomerularFiltrationLPerH = 0.0;  // GFR
    // Enzyme turnover, 1/h, keyed exactly as the pack keys it
    // ("CYP3A4_hepatic", "CYP3A4_intestinal", ...).
    std::map<std::string, double> enzymeDegradation;
    std::string source;   // the pack's own prose citation
    std::string path;     // where it was read from, for the assumption list
    std::string error;    // why it could not be read, when loaded == false
};

// Parsed once and cached; the pack is read-only shipped data, so re-reading it per
// call would be pure I/O for an identical answer.
const Physiology& physiology();

// kdeg for an enzyme, or -1 when the pack does not carry it. Returning -1 rather
// than a plausible 0.02 is deliberate: an unknown turnover makes a time-dependent
// inhibition prediction NotComputed, not approximate.
double enzymeDegradationRate(const std::string& enzymeKey);

}  // namespace biocad::core
