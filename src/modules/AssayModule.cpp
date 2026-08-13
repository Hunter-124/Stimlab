#include "modules/AssayModule.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "assay/Biophysics.h"
#include "assay/Dataset.h"
#include "assay/Design.h"
#include "assay/Fits.h"
#include "assay/Qc.h"

namespace biocad {
namespace {

bool isDoseResponseModel(AssayModel m) {
    switch (m) {
        case AssayModel::FourParameterLogistic:
        case AssayModel::FiveParameterLogistic:
        case AssayModel::MichaelisMenten:
        case AssayModel::Hill:
        case AssayModel::SubstrateInhibition:
        case AssayModel::MorrisonTightBinding:
            return true;
        default:
            return false;
    }
}

// A sensorgram set out of a well list: one curve per analyte concentration, each
// ordered in time.
//
// The injection stop is not a well field, so it is INFERRED as the time of the
// maximum response in the curve, and the inference is reported as a warning on the
// fit rather than presented as metadata. A wrong dissociation start biases kd, so
// the reader has to know the number was inferred.
assay::KineticExperiment kineticFromWells(const std::vector<Well>& wells,
                                          std::vector<std::string>& warnings) {
    assay::KineticExperiment exp;
    if (!wells.empty()) exp.seriesId = wells.front().seriesId;

    std::map<double, std::vector<const Well*>> byConc;
    for (const auto& w : wells) {
        if (w.excluded) continue;
        byConc[w.concentration].push_back(&w);
    }
    for (auto& [conc, pts] : byConc) {
        std::sort(pts.begin(), pts.end(),
                  [](const Well* a, const Well* b) { return a->timeS < b->timeS; });
        assay::KineticCurve c;
        c.concentrationM = conc;
        double peak = -std::numeric_limits<double>::infinity();
        for (const Well* w : pts) {
            c.timeS.push_back(w->timeS);
            c.responseRu.push_back(w->readout);
            if (w->readout > peak) {
                peak = w->readout;
                c.dissociationStartS = w->timeS;
            }
        }
        exp.curves.push_back(std::move(c));
    }
    warnings.push_back("the injection stop is not a well field, so the dissociation start of "
                       "each curve was inferred as the time of its maximum response; supply an "
                       "explicit stop time if the injection did not end at the peak");
    return exp;
}

assay::MeltCurve meltFromWells(const std::vector<Well>& wells) {
    assay::MeltCurve mc;
    if (!wells.empty()) mc.seriesId = wells.front().seriesId;
    std::vector<const Well*> pts;
    for (const auto& w : wells)
        if (!w.excluded) pts.push_back(&w);
    std::sort(pts.begin(), pts.end(),
              [](const Well* a, const Well* b) { return a->temperatureC < b->temperatureC; });
    for (const Well* w : pts) {
        mc.temperatureC.push_back(w->temperatureC);
        mc.signal.push_back(w->readout);
    }
    return mc;
}

FitResult refuse(AssayModel model, const std::string& seriesId, const std::string& missing) {
    FitResult f;
    f.model = model;
    f.seriesId = seriesId;
    f.converged = false;
    f.note = "not fitted: " + missing;
    f.warnings.push_back(f.note);
    return f;
}

}  // namespace

std::optional<AssayDataset> RealAssay::import(const std::string& text) const {
    std::string err;
    AssayDataset ds = assay::importText(text, &err);
    // The engine reports a hard failure through `err` and a recoverable one through
    // AssayDataset::warnings. std::nullopt is reserved for "this was not tabular at
    // all", so a plate that imported with complaints still reaches the panel.
    if (ds.plates.empty()) {
        if (!err.empty()) return std::nullopt;
        return std::nullopt;
    }
    return ds;
}

QcReport RealAssay::qc(const Plate& p) const { return assay::plateQc(p); }

FitResult RealAssay::fit(const std::vector<Well>& series, AssayModel model, bool robust) const {
    assay::FitOptions opts;
    opts.robust = robust;

    if (isDoseResponseModel(model)) return assay::fitSeries(series, model, opts);

    const std::string sid = series.empty() ? std::string() : series.front().seriesId;
    switch (model) {
        case AssayModel::LangmuirKinetics:
        case AssayModel::MassTransportKinetics: {
            std::vector<std::string> warnings;
            const assay::KineticExperiment exp = kineticFromWells(series, warnings);
            if (exp.curves.empty())
                return refuse(model, sid, "no sensorgram points survived exclusion");
            FitResult f = model == AssayModel::LangmuirKinetics
                              ? assay::fitLangmuirKinetics(exp)
                              : assay::fitMassTransportKinetics(exp);
            for (const auto& w : warnings) f.warnings.push_back(w);
            return f;
        }
        case AssayModel::BoltzmannMelt:
        case AssayModel::TwoStateThermodynamic: {
            const assay::MeltCurve mc = meltFromWells(series);
            if (mc.temperatureC.size() < 5)
                return refuse(model, sid,
                              "a melt needs at least five temperature points with a "
                              "temperature_c column");
            return model == AssayModel::BoltzmannMelt ? assay::fitBoltzmannMelt(mc)
                                                      : assay::fitTwoStateMelt(mc);
        }
        case AssayModel::WisemanIsotherm:
            // Deliberate refusal, not an omission: n, K and dH are only meaningful
            // against a known cell volume, macromolecule and titrant concentration
            // and a blank titration, and none of those is a property of a well.
            return refuse(model, sid,
                          "an ITC isotherm needs the cell volume, the macromolecule and "
                          "titrant concentrations and the blank heat of dilution, which are "
                          "experiment metadata and not well fields; load an ITC experiment "
                          "rather than a well list");
        default:
            return refuse(model, sid, "no fitter is registered for that model");
    }
}

ModelComparison RealAssay::compare(const std::vector<Well>&       series,
                                   const std::vector<AssayModel>& candidates) const {
    // Only the series-fittable candidates are forwarded: ranking a sensorgram model
    // against a 4PL by AICc over a concentration series would compare fits to
    // different data.
    std::vector<AssayModel> usable;
    for (AssayModel m : candidates)
        if (isDoseResponseModel(m)) usable.push_back(m);
    ModelComparison c = assay::compareModels(series, usable, {});
    if (usable.size() != candidates.size())
        c.conclusion += " (trace-based models were dropped from the ranking: they are not "
                        "fitted to a concentration series)";
    return c;
}

ModelComparison RealAssay::inhibitionModality(const std::vector<Well>& matrix) const {
    return assay::fitInhibitionModalityFromWells(matrix, {});
}

AssayDesignReport RealAssay::simulate(const AssayDesignSpec& spec) const {
    return assay::simulateDesign(spec);
}

}  // namespace biocad
