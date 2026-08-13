// chem/Solubility.h - buffer capacity, pH-dependent solubility, the BCS numbers,
// and the dissolution/precipitation time course.
//
// WHY THIS FILE IS SHAPED THE WAY IT IS. Every one of these calculations is exact
// physics on top of numbers that must be MEASURED: an intrinsic solubility needs a
// melting point (or a measured S0), a salt plateau needs a Ksp and a counterion
// concentration, a precipitation rate needs a rate constant. None of those are
// predictable from a structure, so none of them are defaulted. When one is absent
// the dependent Quantity is notComputed() naming it, and the curve simply does not
// have that branch - drawing a salt kink whose position is unknown would be
// fabricating the single most important feature of the plot.
//
// VAN SLYKE. beta = 2.303 (Kw/[H+] + [H+] + sum_i C_i Ka_i [H+] / (Ka_i + [H+])^2)
// per Van Slyke's definition of buffer value, dC_base/dpH. The water terms are
// never dropped: they are what makes unbuffered water's beta at pH 7 nonzero, and
// a "buffer capacity" that reads exactly zero for pure water is wrong.
//
// GENERAL SOLUBILITY EQUATION. log S0 = 0.5 - 0.01 (MP_C - 25) - logP, the Jain &
// Yalkowsky form (their revision of the older Yalkowsky-Valvani expression, which
// carried an entropy-of-melting term). It is Provenance::Predicted and therefore
// must carry an error bar; the published average absolute error is about 0.4-0.5
// log10 units on drug-like sets, which is what gseCitation() states and what the
// returned Quantity's error bar is propagated from. No DOI is asserted here
// because the exact identifier was not verified in this environment.
//
// SPECIATION. Only the monoprotic acid/base cases are needed for a pH-solubility
// profile, and they are closed-form Henderson-Hasselbalch, so they are implemented
// locally rather than routed through the general component-tableau solver in
// chem/Speciation.* (which did not yet exist in the tree when this file was
// written). Nothing here duplicates that solver's polyprotic capability.
#pragma once

#include <string>
#include <vector>

#include "data/Ionization.h"

namespace biocad::chem {

// --------------------------------------------------------------------- buffers

struct BufferSpec {
    std::vector<BufferComponent> components;   // conjugate pairs, pKa + total molarity
    double pHMin = 1.0;
    double pHMax = 13.0;
    double pHStep = 0.02;
    double kw = 1.0e-14;   // water ion product at the stated temperature
};

// Van Slyke buffer value at one pH, mol/L per pH unit. Exposed because it is the
// phase's headline fixture and a test should be able to hit it without a curve.
double vanSlykeBeta(const std::vector<BufferComponent>& components, double pH,
                    double kw = 1.0e-14);

// The curve plus beta at pH 7.4 and the maximum with the pH it occurs at. The
// maximum is located on the sampled grid, so its pH is only as precise as pHStep;
// that limitation is recorded in BufferReport::assumptions with the real number.
BufferReport bufferCapacity(const BufferSpec& spec);

// ------------------------------------------------------------------ solubility

enum class IonizationKind {
    Neutral,          // no ionizable group: S(pH) is flat
    MonoproticAcid,   // HA <-> A- + H+
    MonoproticBase    // BH+ <-> B + H+
};

// Everything the pH-solubility profile can use. Each `has*` flag exists so that
// "absent" is distinguishable from "zero" - a Ksp of 0 is a real (insoluble-salt)
// statement, whereas an absent Ksp means there is no salt branch to draw.
struct SolubilityInput {
    std::string    moleculeId;
    IonizationKind kind = IonizationKind::Neutral;

    double pKa = 0.0;              bool hasPKa = false;
    double logP = 0.0;             bool hasLogP = false;
    double meltingPointC = 0.0;    bool hasMeltingPoint = false;

    // A measured S0 always wins over the GSE estimate and carries its own source.
    double      measuredS0Molar = 0.0;  bool hasMeasuredS0 = false;
    std::string measuredS0Source;

    // Salt branch. `ksp` is the 1:1 salt solubility product in (mol/L)^2 and
    // `counterionMolar` is the common-ion concentration that suppresses it.
    double ksp = 0.0;              bool hasKsp = false;
    double counterionMolar = 0.0;  bool hasCounterion = false;

    // BCS inputs.
    double doseMg = 0.0;                bool hasDose = false;
    double molWeight = 0.0;             bool hasMolWeight = false;
    double particleRadiusUm = 0.0;      bool hasParticleRadius = false;
    double diffusivityCm2PerS = 0.0;    bool hasDiffusivity = false;
    double densityGPerCm3 = 0.0;        bool hasDensity = false;
    double peffCmPerS = 0.0;            bool hasPeff = false;

    double pHMin = 1.0;
    double pHMax = 10.0;
    double pHStep = 0.02;
};

// Intrinsic solubility S0 in mol/L from the General Solubility Equation. The
// caller must have a melting point AND a logP; this overload assumes both are
// present, and phSolubility() is what decides whether to call it.
Quantity gseIntrinsicSolubility(double logP, double meltingPointC);

// The exact source string every GSE-derived Quantity carries, error bar included.
const char* gseCitation();

// pH-solubility with the salt-limited plateau, plus Do/Dn/An. Curve points are
// flagged saltLimited on the plateau side of pHmax, so the renderer never has to
// re-derive where the kink was.
SolubilityReport phSolubility(const SolubilityInput& in);

// ----------------------------------------------------------------- dissolution

struct DissolutionInput {
    double doseMg = 0.0;
    double molWeight = 0.0;             // g/mol, to convert mg <-> mol/L
    double volumeL = 0.250;             // BCS reference gastric volume
    double initialRadiusUm = 0.0;
    double densityGPerCm3 = 0.0;
    double diffusivityCm2PerS = 0.0;
    double diffusionLayerUm = 0.0;      // Noyes-Whitney film thickness h
    double solubilityMolar = 0.0;       // Cs seen by the dissolving surface

    // pH-shift precipitation: dC/dt = -kppt (C - S). kppt is supplied or fitted,
    // never predicted, so `hasKppt` false means no precipitation term at all.
    bool   precipitation = false;
    double kpptPerS = 0.0;              bool hasKppt = false;
    double precipSolubilityMolar = 0.0; // S the solution relaxes back down to

    double horizonS = 3600.0;
    double stepS = 0.5;
};

// Noyes-Whitney dissolution of Hixson-Crowell shrinking spheres, integrated with
// numeric::rk4Integrate. Solid, dissolved and precipitated mass are three separate
// states whose derivatives sum to zero identically, so their total is conserved to
// floating-point roundoff; DissolutionReport::maxMassImbalance reports the worst
// observed deviation in mg rather than asking to be trusted.
DissolutionReport dissolutionTimeCourse(const DissolutionInput& in);

// Least-squares fit of kppt to an observed concentration decay (times in s,
// concentrations in mol/L) against dC/dt = -kppt (C - S). Returns notComputed when
// there is nothing to fit. This is the ONLY sanctioned origin of a kppt besides
// the user typing one.
Quantity fitPrecipitationRate(const std::vector<double>& timeS,
                              const std::vector<double>& concentrationMolar,
                              double solubilityMolar);

}  // namespace biocad::chem
