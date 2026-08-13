// sim/Sbml.h - SBML Core reader and writer for the rate laws BioCAD implements.
//
// WHY NOT libSBML. libSBML is LGPL. BioCAD ships as a single fully static Windows
// binary (the windows-static preset links the static CRT and every dependency), and
// the LGPL's relinking obligation is incompatible with distributing that binary
// without also shipping the object files needed to relink it. So libSBML is NOT
// linked, statically or otherwise. The subset BioCAD needs is a few hundred lines of
// pugixml (MIT) over a well-specified document, and that is what this file is.
//
// WHAT IS SUPPORTED. SBML Core: model, listOfCompartments, listOfSpecies,
// listOfParameters, listOfReactions with listOfReactants / listOfProducts /
// listOfModifiers, and a kineticLaw whose Content MathML is one of the five rate
// laws sim::Network can integrate (mass action, reversible mass action,
// Michaelis-Menten, reversible Michaelis-Menten, Hill). Level 3 Version 2 is the
// export target; Level 2 is also READ, because the curated CC0 BioModels corpus is
// overwhelmingly Level 2 and a Level-3-only reader could not open a single committed
// fixture.
//
// WHAT IS REFUSED, BY NAME. Function definitions, assignment/rate/algebraic rules,
// events, constraints, initial assignments, and any kinetic law outside those five
// forms are reported through the `error` out-parameter naming the construct and the
// element it appeared on. Silently dropping an assignment rule would change the
// model's meaning and still produce a curve, which is the worst possible outcome:
// importSbml returns std::nullopt instead.
#pragma once

#include <optional>
#include <string>

#include "data/Systems.h"

namespace biocad::sim {

// Parses an SBML document. Returns std::nullopt with `error` naming the first
// unsupported construct (element type plus the id it appeared on).
std::optional<NetworkSpec> importSbml(const std::string& xml, std::string* error);

// Writes SBML Level 3 Version 2 Core. Every rate law is emitted in the canonical
// Content MathML form importSbml recognises, so export -> import is structurally
// lossless for anything BioCAD can represent.
std::string exportSbml(const NetworkSpec& network);

}  // namespace biocad::sim
