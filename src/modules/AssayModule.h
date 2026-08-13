// modules/AssayModule.h - the Phase 10 adapter between the assay engines and the
// rest of BioCAD.
//
// WHAT THIS FILE IS FOR. src/assay/{Dataset,Qc,Fits,Biophysics,Design}.* know how
// to import a plate, judge it, fit it and simulate it, but they know nothing about
// panels, tools or services. This adapter is the single dispatch point that decides
// WHICH engine a requested AssayModel belongs to: the dose-response and enzyme
// models are fitted from a concentration series, and the SPR/DSF/ITC models are
// fitted from a trace. That distinction is not cosmetic - a sensorgram is a
// response-versus-time curve per analyte concentration, and a melt is signal versus
// temperature, so a std::vector<Well> has to be reshaped into the experiment type
// each engine actually takes.
//
// WHAT IT REFUSES. An ITC isotherm needs the cell volume, the macromolecule and
// titrant concentrations and the blank titration; none of those are properties of a
// well, so an ITC fit requested from a bare well list comes back converged = false
// naming what was missing. That is the honest answer, and inventing a cell volume
// to make the call succeed would produce a thermodynamic parameter set from a
// number nobody measured.
#pragma once

#include <string>
#include <vector>

#include "contracts/IModules.h"

namespace biocad {

class RealAssay final : public IAssayModule {
public:
    std::optional<AssayDataset> import(const std::string& text) const override;
    QcReport                    qc(const Plate& p) const override;
    FitResult                   fit(const std::vector<Well>& series, AssayModel model,
                                    bool robust) const override;
    ModelComparison             compare(const std::vector<Well>&            series,
                                        const std::vector<AssayModel>&      candidates) const override;
    ModelComparison             inhibitionModality(const std::vector<Well>& matrix) const override;
    AssayDesignReport           simulate(const AssayDesignSpec& spec) const override;
};

}  // namespace biocad
