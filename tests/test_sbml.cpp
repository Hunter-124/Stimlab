// tests/test_sbml.cpp - the SBML Core subset, against three REAL committed CC0
// BioModels documents and against BioCAD's own writer.
//
// The fixtures are the genuine article, not hand-written stand-ins:
//   BIOMD0000000052  Brands 2002, monosaccharide-casein Maillard kinetics (11 species)
//   BIOMD0000000050  Martins 2003, Amadori degradation (14 species)
//   BIOMD0000000035  Vilar 2002, circadian oscillator (10 species)
// All three are curated BioModels entries, released CC0, and all three are SBML
// Level 2 - which is why the reader accepts level 2 as well as level 3. Nothing in
// this file touches the network.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "sim/Sbml.h"
#include "sim/Solvers.h"

using namespace biocad;
using namespace biocad::sim;

namespace {

std::string fixture(const std::string& name) {
    std::ifstream in(std::string(BIOCAD_TEST_FIXTURES) + "/sbml/" + name, std::ios::binary);
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

}  // namespace

TEST_CASE("the three committed CC0 BioModels fixtures import and integrate", "[sbml]") {
    struct Expect { const char* file; std::size_t species; std::size_t reactions; };
    const Expect cases[] = {{"BIOMD0000000052.xml", 11, 11},
                            {"BIOMD0000000050.xml", 14, 16},
                            {"BIOMD0000000035.xml", 10, 16}};
    for (const Expect& c : cases) {
        const std::string xml = fixture(c.file);
        REQUIRE_FALSE(xml.empty());
        std::string error;
        const auto net = importSbml(xml, &error);
        REQUIRE(net.has_value());
        REQUIRE(error.empty());
        REQUIRE(net->species.size() == c.species);
        REQUIRE(net->reactions.size() == c.reactions);
        // Importing runs the structural analysis, so the conserved moieties of a real
        // published model are available immediately.
        REQUIRE(net->conservationLaws.size() == net->conservationLabels.size());

        IntegrationOptions o;
        o.horizon = 1.0;
        o.outputPoints = 11;
        o.relativeTolerance = 1e-8;
        o.absoluteTolerance = 1e-14;
        const TimeCourse tc = integrate(*net, o);
        REQUIRE(tc.times.size() == 11);
        REQUIRE(tc.trajectories.size() == c.species);
        REQUIRE(tc.worstConservationDrift < 1e-10);
    }
}

TEST_CASE("a real fixture round-trips structurally through the writer", "[sbml]") {
    std::string error;
    const auto original = importSbml(fixture("BIOMD0000000052.xml"), &error);
    REQUIRE(original.has_value());
    std::string secondError;
    const auto reread = importSbml(exportSbml(*original), &secondError);
    REQUIRE(reread.has_value());
    REQUIRE(original->species.size() == reread->species.size());
    REQUIRE(original->reactions.size() == reread->reactions.size());
    for (std::size_t i = 0; i < original->species.size(); ++i) {
        REQUIRE(original->species[i].id == reread->species[i].id);
        REQUIRE(original->species[i].initialConcentration ==
                reread->species[i].initialConcentration);
        REQUIRE(original->species[i].boundary == reread->species[i].boundary);
    }
    for (std::size_t i = 0; i < original->reactions.size(); ++i) {
        REQUIRE(original->reactions[i].id == reread->reactions[i].id);
        REQUIRE(original->reactions[i].law == reread->reactions[i].law);
        REQUIRE(original->reactions[i].reactants == reread->reactions[i].reactants);
        REQUIRE(original->reactions[i].products == reread->reactions[i].products);
        REQUIRE(original->reactions[i].parameters == reread->reactions[i].parameters);
    }
    // Structural identity is not enough on its own: the two must also integrate to
    // the same numbers, which is what proves no rate constant was quietly rescaled.
    IntegrationOptions o;
    o.horizon = 5.0;
    o.outputPoints = 6;
    o.relativeTolerance = 1e-9;
    o.absoluteTolerance = 1e-14;
    const TimeCourse a = integrate(*original, o), b = integrate(*reread, o);
    for (std::size_t i = 0; i < a.trajectories.size(); ++i)
        for (std::size_t k = 0; k < a.times.size(); ++k)
            REQUIRE(a.trajectories[i][k] == b.trajectories[i][k]);
}

TEST_CASE("every unsupported SBML construct is refused BY NAME", "[sbml]") {
    struct Case { const char* xml; const char* mustName; };
    const Case cases[] = {
        {R"(<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version2/core"
            level="3" version="2"><model id="m">
            <listOfSpecies><species id="A" compartment="c" initialConcentration="1"/></listOfSpecies>
            <listOfRules><assignmentRule variable="A"><math/></assignmentRule></listOfRules>
            <listOfReactions/></model></sbml>)", "assignmentRule"},
        {R"(<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version2/core"
            level="3" version="2"><model id="m">
            <listOfEvents><event id="fire"/></listOfEvents></model></sbml>)", "event"},
        {R"(<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version2/core"
            level="3" version="2"><model id="m">
            <listOfFunctionDefinitions><functionDefinition id="f"/></listOfFunctionDefinitions>
            </model></sbml>)", "functionDefinition"},
        {R"(<?xml version="1.0"?><sbml xmlns="x" level="1" version="2"><model id="m"/></sbml>)",
         "level 1"},
        // A well-formed law that uses exp(): named by OPERATOR, never swallowed.
        {R"(<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version2/core"
            level="3" version="2"><model id="m">
            <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
            <listOfSpecies><species id="A" compartment="c" initialConcentration="1"/>
            <species id="B" compartment="c" initialConcentration="0"/></listOfSpecies>
            <listOfReactions><reaction id="R1" reversible="false">
            <listOfReactants><speciesReference species="A" stoichiometry="1"/></listOfReactants>
            <listOfProducts><speciesReference species="B" stoichiometry="1"/></listOfProducts>
            <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><exp/><ci>A</ci></apply></math>
            <listOfLocalParameters><localParameter id="k" value="1"/></listOfLocalParameters>
            </kineticLaw></reaction></listOfReactions></model></sbml>)", "exp"},
        // Mass-action SHAPED, but over a species that is not a reactant - the
        // enzyme-modified form real reconstructions use. Refused, not approximated.
        {R"(<?xml version="1.0"?><sbml xmlns="http://www.sbml.org/sbml/level3/version2/core"
            level="3" version="2"><model id="m">
            <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
            <listOfSpecies><species id="A" compartment="c" initialConcentration="1"/>
            <species id="B" compartment="c" initialConcentration="0"/>
            <species id="E" compartment="c" initialConcentration="2"/></listOfSpecies>
            <listOfReactions><reaction id="Rbad" reversible="false">
            <listOfReactants><speciesReference species="A" stoichiometry="1"/></listOfReactants>
            <listOfProducts><speciesReference species="B" stoichiometry="1"/></listOfProducts>
            <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>A</ci><ci>E</ci></apply></math>
            <listOfLocalParameters><localParameter id="k" value="1"/></listOfLocalParameters>
            </kineticLaw></reaction></listOfReactions></model></sbml>)", "Rbad"},
    };
    for (const Case& c : cases) {
        std::string error;
        const auto net = importSbml(c.xml, &error);
        REQUIRE_FALSE(net.has_value());
        REQUIRE(error.find(c.mustName) != std::string::npos);
    }
}

TEST_CASE("all five rate laws survive a write and read cycle", "[sbml]") {
    NetworkSpec n;
    n.id = "all-laws";
    for (const char* id : {"A", "B", "C", "D", "E", "F"}) {
        SpeciesSpec s;
        s.id = id;
        s.name = id;
        s.compartment = "cell";
        s.initialConcentration = 1.0;
        n.species.push_back(s);
    }
    auto rx = [](const char* id, RateLaw law, std::vector<double> parameters,
                 std::vector<std::pair<std::string, double>> reactants,
                 std::vector<std::pair<std::string, double>> products) {
        ReactionSpec r;
        r.id = id;
        r.law = law;
        r.parameters = std::move(parameters);
        r.reactants = std::move(reactants);
        r.products = std::move(products);
        r.reversible = law == RateLaw::ReversibleMassAction ||
                       law == RateLaw::ReversibleMichaelisMenten;
        return r;
    };
    n.reactions = {
        rx("ma", RateLaw::MassAction, {0.31}, {{"A", 2}}, {{"B", 1}}),
        rx("rma", RateLaw::ReversibleMassAction, {0.7, 0.25}, {{"B", 1}}, {{"C", 1}}),
        rx("mm", RateLaw::MichaelisMenten, {2.5, 0.4}, {{"C", 1}}, {{"D", 1}}),
        rx("hill", RateLaw::Hill, {1.75, 0.6, 2.0}, {{"D", 1}}, {{"E", 1}}),
        rx("rmm", RateLaw::ReversibleMichaelisMenten, {1.2, 0.3, 0.8, 0.5}, {{"E", 1}},
           {{"F", 1}}),
    };
    std::string error;
    const auto back = importSbml(exportSbml(n), &error);
    REQUIRE(back.has_value());
    REQUIRE(back->reactions.size() == 5);
    for (std::size_t i = 0; i < 5; ++i) {
        REQUIRE(n.reactions[i].law == back->reactions[i].law);
        REQUIRE(n.reactions[i].parameters == back->reactions[i].parameters);
        REQUIRE(n.reactions[i].reactants == back->reactions[i].reactants);
        REQUIRE(n.reactions[i].products == back->reactions[i].products);
    }
}
