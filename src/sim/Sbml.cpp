#include "sim/Sbml.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <pugixml.hpp>

#include "sim/Network.h"

namespace biocad::sim {
namespace {

// SBML documents normally use a default namespace, so pugixml sees unprefixed names -
// but some producers write sbml:model. Matching on the LOCAL name handles both without
// a namespace-aware parser.
std::string_view localName(const char* name) {
    std::string_view s(name);
    const std::size_t colon = s.rfind(':');
    return colon == std::string_view::npos ? s : s.substr(colon + 1);
}

pugi::xml_node child(const pugi::xml_node& n, std::string_view name) {
    for (pugi::xml_node c : n.children())
        if (localName(c.name()) == name) return c;
    return {};
}

std::vector<pugi::xml_node> children(const pugi::xml_node& n, std::string_view name) {
    std::vector<pugi::xml_node> out;
    for (pugi::xml_node c : n.children())
        if (localName(c.name()) == name) out.push_back(c);
    return out;
}

std::string trim(std::string s) {
    const char* ws = " \t\r\n";
    const std::size_t a = s.find_first_not_of(ws);
    if (a == std::string::npos) return {};
    const std::size_t b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

std::string fmt(double v) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}

// ------------------------------------------------------------ MathML matching
// A tiny view over the Content MathML subset the five rate laws need. Everything
// else is refused by name, so this never has to be a general expression evaluator.
struct Symbols {
    // id -> value for parameters (local shadowing global), and the set of species ids.
    std::unordered_map<std::string, double> parameters;
    std::unordered_map<std::string, std::size_t> speciesIndex;
    std::unordered_map<std::string, double> compartmentSize;
};

enum class TokenKind { Unknown, Number, Parameter, Species, Compartment };

struct Token {
    TokenKind kind = TokenKind::Unknown;
    double    number = 0;
    std::string id;
};

Token classify(const pugi::xml_node& n, const Symbols& sym) {
    Token t;
    const std::string_view name = localName(n.name());
    if (name == "cn") {
        t.kind = TokenKind::Number;
        t.number = std::atof(trim(n.text().as_string()).c_str());
        return t;
    }
    if (name == "ci") {
        t.id = trim(n.text().as_string());
        if (sym.speciesIndex.count(t.id)) t.kind = TokenKind::Species;
        else if (sym.parameters.count(t.id)) t.kind = TokenKind::Parameter;
        else if (sym.compartmentSize.count(t.id)) t.kind = TokenKind::Compartment;
        return t;
    }
    return t;
}

// The operator of an <apply>, e.g. "times"; empty for a non-apply node.
std::string_view applyOperator(const pugi::xml_node& n) {
    if (localName(n.name()) != "apply") return {};
    for (pugi::xml_node c : n.children())
        if (c.type() == pugi::node_element) return localName(c.name());
    return {};
}

std::vector<pugi::xml_node> applyArguments(const pugi::xml_node& n) {
    std::vector<pugi::xml_node> out;
    bool first = true;
    for (pugi::xml_node c : n.children()) {
        if (c.type() != pugi::node_element) continue;
        if (first) { first = false; continue; }   // the operator itself
        out.push_back(c);
    }
    return out;
}

// A product of one rate constant and species powers: k * prod(S_i^a_i).
struct Product {
    bool                                 ok = false;
    double                               constant = 1.0;
    std::string                          constantId;
    std::map<std::string, double>        powers;      // species id -> exponent
};

Product matchProduct(const pugi::xml_node& n, const Symbols& sym) {
    Product p;
    std::vector<pugi::xml_node> factors;
    if (applyOperator(n) == "times") factors = applyArguments(n);
    else factors.push_back(n);   // a bare k or a bare species is a one-factor product

    bool sawConstant = false;
    for (const pugi::xml_node& f : factors) {
        const std::string_view op = applyOperator(f);
        if (op == "power") {
            const auto args = applyArguments(f);
            if (args.size() != 2) return p;
            const Token base = classify(args[0], sym);
            const Token exp = classify(args[1], sym);
            if (base.kind != TokenKind::Species) return p;
            double e = 0;
            if (exp.kind == TokenKind::Number) e = exp.number;
            else if (exp.kind == TokenKind::Parameter) e = sym.parameters.at(exp.id);
            else return p;
            p.powers[base.id] += e;
            continue;
        }
        if (!op.empty()) return p;   // any other nested operator: not a plain product
        const Token t = classify(f, sym);
        switch (t.kind) {
            case TokenKind::Species: p.powers[t.id] += 1.0; break;
            case TokenKind::Parameter:
                if (sawConstant) return p;   // two rate constants is not this law
                p.constant *= sym.parameters.at(t.id);
                p.constantId = t.id;
                sawConstant = true;
                break;
            case TokenKind::Number: p.constant *= t.number; break;
            case TokenKind::Compartment:
                // A compartment id as a multiplicative factor converts a
                // concentration rate to an amount rate. It is a no-op only when the
                // compartment has unit size; a non-unit size is refused upstream.
                p.constant *= sym.compartmentSize.at(t.id);
                break;
            case TokenKind::Unknown: return p;
        }
    }
    p.ok = sawConstant || !p.powers.empty();
    return p;
}

bool powersMatch(const std::map<std::string, double>& powers,
                 const std::vector<std::pair<std::string, double>>& side) {
    std::map<std::string, double> want;
    for (const auto& [id, s] : side) want[id] += s;
    if (want.size() != powers.size()) return false;
    for (const auto& [id, s] : want) {
        auto it = powers.find(id);
        if (it == powers.end() || std::abs(it->second - s) > 1e-9) return false;
    }
    return true;
}

// divide(times(V, S), plus(K, S))  ->  Michaelis-Menten in S.
bool matchMichaelisMenten(const pugi::xml_node& n, const Symbols& sym, const std::string& substrate,
                          double& vmax, double& km) {
    if (applyOperator(n) != "divide") return false;
    const auto args = applyArguments(n);
    if (args.size() != 2) return false;
    const Product num = matchProduct(args[0], sym);
    if (!num.ok || num.powers.size() != 1) return false;
    if (num.powers.begin()->first != substrate) return false;
    if (std::abs(num.powers.begin()->second - 1.0) > 1e-12) return false;
    if (applyOperator(args[1]) != "plus") return false;
    const auto den = applyArguments(args[1]);
    if (den.size() != 2) return false;
    double kmValue = 0;
    bool sawSubstrate = false, sawKm = false;
    for (const pugi::xml_node& d : den) {
        const Token t = classify(d, sym);
        if (t.kind == TokenKind::Species && t.id == substrate) sawSubstrate = true;
        else if (t.kind == TokenKind::Parameter) { kmValue = sym.parameters.at(t.id); sawKm = true; }
        else if (t.kind == TokenKind::Number) { kmValue = t.number; sawKm = true; }
        else return false;
    }
    if (!sawSubstrate || !sawKm) return false;
    vmax = num.constant;
    km = kmValue;
    return true;
}

// divide(times(V, power(S,n)), plus(power(K,n), power(S,n)))  ->  Hill.
bool matchHill(const pugi::xml_node& n, const Symbols& sym, const std::string& substrate,
               double& vmax, double& k, double& hill) {
    if (applyOperator(n) != "divide") return false;
    const auto args = applyArguments(n);
    if (args.size() != 2) return false;
    const Product num = matchProduct(args[0], sym);
    if (!num.ok || num.powers.size() != 1 || num.powers.begin()->first != substrate) return false;
    const double exponent = num.powers.begin()->second;
    if (applyOperator(args[1]) != "plus") return false;
    const auto den = applyArguments(args[1]);
    if (den.size() != 2) return false;
    bool sawSubstratePower = false, sawKPower = false;
    for (const pugi::xml_node& d : den) {
        if (applyOperator(d) != "power") return false;
        const auto pa = applyArguments(d);
        if (pa.size() != 2) return false;
        const Token base = classify(pa[0], sym);
        const Token exp = classify(pa[1], sym);
        double e = 0;
        if (exp.kind == TokenKind::Number) e = exp.number;
        else if (exp.kind == TokenKind::Parameter) e = sym.parameters.at(exp.id);
        else return false;
        if (std::abs(e - exponent) > 1e-9) return false;
        if (base.kind == TokenKind::Species && base.id == substrate) sawSubstratePower = true;
        else if (base.kind == TokenKind::Parameter) { k = sym.parameters.at(base.id); sawKPower = true; }
        else if (base.kind == TokenKind::Number) { k = base.number; sawKPower = true; }
        else return false;
    }
    if (!sawSubstratePower || !sawKPower) return false;
    vmax = num.constant;
    hill = exponent;
    return true;
}

// divide(x, y) with x = ci/param -> a plain ratio; used by the reversible MM form.
bool matchRatio(const pugi::xml_node& n, const Symbols& sym, std::string& speciesId, double& km) {
    if (applyOperator(n) != "divide") return false;
    const auto args = applyArguments(n);
    if (args.size() != 2) return false;
    const Token s = classify(args[0], sym);
    const Token k = classify(args[1], sym);
    if (s.kind != TokenKind::Species) return false;
    if (k.kind == TokenKind::Parameter) km = sym.parameters.at(k.id);
    else if (k.kind == TokenKind::Number) km = k.number;
    else return false;
    speciesId = s.id;
    return true;
}

// divide(minus(divide(times(Vf,S),Ks), divide(times(Vr,P),Kp)),
//        plus(1, divide(S,Ks), divide(P,Kp)))
bool matchReversibleMm(const pugi::xml_node& n, const Symbols& sym, const std::string& substrate,
                       const std::string& product, double& vf, double& ks, double& vr,
                       double& kp) {
    if (applyOperator(n) != "divide") return false;
    const auto args = applyArguments(n);
    if (args.size() != 2) return false;
    if (applyOperator(args[0]) != "minus") return false;
    const auto num = applyArguments(args[0]);
    if (num.size() != 2) return false;
    auto limb = [&](const pugi::xml_node& x, const std::string& want, double& v, double& k) {
        if (applyOperator(x) != "divide") return false;
        const auto a = applyArguments(x);
        if (a.size() != 2) return false;
        const Product p = matchProduct(a[0], sym);
        if (!p.ok || p.powers.size() != 1 || p.powers.begin()->first != want) return false;
        const Token kt = classify(a[1], sym);
        if (kt.kind == TokenKind::Parameter) k = sym.parameters.at(kt.id);
        else if (kt.kind == TokenKind::Number) k = kt.number;
        else return false;
        v = p.constant;
        return true;
    };
    if (!limb(num[0], substrate, vf, ks)) return false;
    if (!limb(num[1], product, vr, kp)) return false;
    if (applyOperator(args[1]) != "plus") return false;
    const auto den = applyArguments(args[1]);
    if (den.size() != 3) return false;
    bool sawOne = false, sawS = false, sawP = false;
    for (const pugi::xml_node& d : den) {
        const Token t = classify(d, sym);
        if (t.kind == TokenKind::Number && std::abs(t.number - 1.0) < 1e-12) { sawOne = true; continue; }
        std::string id;
        double k = 0;
        if (!matchRatio(d, sym, id, k)) return false;
        if (id == substrate && std::abs(k - ks) < 1e-9) sawS = true;
        else if (id == product && std::abs(k - kp) < 1e-9) sawP = true;
        else return false;
    }
    return sawOne && sawS && sawP;
}

}  // namespace

std::optional<NetworkSpec> importSbml(const std::string& xml, std::string* error) {
    auto fail = [&](const std::string& why) -> std::optional<NetworkSpec> {
        if (error) *error = why;
        return std::nullopt;
    };
    pugi::xml_document doc;
    const pugi::xml_parse_result parsed = doc.load_buffer(xml.data(), xml.size());
    if (!parsed) return fail(std::string("XML is not well formed: ") + parsed.description());

    pugi::xml_node sbml = child(doc, "sbml");
    if (!sbml) return fail("document has no <sbml> root element");
    const int level = sbml.attribute("level").as_int(0);
    const int version = sbml.attribute("version").as_int(0);
    if (level != 2 && level != 3)
        return fail("unsupported SBML level " + std::to_string(level) + " version " +
                    std::to_string(version) + " (BioCAD reads Core level 2 and level 3)");

    pugi::xml_node model = child(sbml, "model");
    if (!model) return fail("<sbml> contains no <model>");

    NetworkSpec out;
    out.id = model.attribute("id").as_string(model.attribute("name").as_string("sbml-model"));

    // Constructs that change the model's meaning and have no representation in a
    // sim::Network. Each is refused by name, with the id of the first offender: a
    // dropped assignment rule still produces a plausible curve, which is worse than
    // an error.
    struct Refusal { const char* list; const char* item; const char* what; };
    static const Refusal kRefusals[] = {
        {"listOfFunctionDefinitions", "functionDefinition", "a user-defined function"},
        {"listOfRules", "assignmentRule", "an assignment rule"},
        {"listOfRules", "rateRule", "a rate rule"},
        {"listOfRules", "algebraicRule", "an algebraic rule"},
        {"listOfEvents", "event", "a discrete event"},
        {"listOfConstraints", "constraint", "a constraint"},
        {"listOfInitialAssignments", "initialAssignment", "an initial assignment"},
    };
    for (const Refusal& r : kRefusals) {
        pugi::xml_node list = child(model, r.list);
        if (!list) continue;
        for (const pugi::xml_node& item : children(list, r.item)) {
            std::string who = item.attribute("id").as_string("");
            if (who.empty()) who = item.attribute("variable").as_string("(unnamed)");
            return fail(std::string("unsupported SBML construct <") + r.item + "> (" + r.what +
                        ") on '" + who + "': BioCAD's reaction networks are stoichiometry plus a "
                        "rate law, and importing this document without it would silently change "
                        "the model");
        }
    }

    Symbols sym;
    for (const pugi::xml_node& c : children(child(model, "listOfCompartments"), "compartment")) {
        const std::string id = c.attribute("id").as_string("");
        const double size = c.attribute("size") ? c.attribute("size").as_double()
                                                : c.attribute("volume").as_double(1.0);
        sym.compartmentSize[id] = size == 0 ? 1.0 : size;
        if (std::abs(sym.compartmentSize[id] - 1.0) > 1e-12)
            out.warnings.push_back("compartment '" + id + "' has size " +
                                   fmt(sym.compartmentSize[id]) +
                                   ": BioCAD integrates concentrations, so any compartment factor "
                                   "in a rate law is folded into the rate constant");
    }

    for (const pugi::xml_node& s : children(child(model, "listOfSpecies"), "species")) {
        SpeciesSpec spec;
        spec.id = s.attribute("id").as_string("");
        spec.name = s.attribute("name").as_string(spec.id.c_str());
        spec.compartment = s.attribute("compartment").as_string("");
        if (s.attribute("initialConcentration")) {
            spec.initialConcentration = s.attribute("initialConcentration").as_double();
        } else if (s.attribute("initialAmount")) {
            const double v = sym.compartmentSize.count(spec.compartment)
                                 ? sym.compartmentSize[spec.compartment]
                                 : 1.0;
            spec.initialConcentration = s.attribute("initialAmount").as_double() / v;
        }
        spec.boundary = s.attribute("boundaryCondition").as_bool(false) ||
                        s.attribute("constant").as_bool(false);
        sym.speciesIndex.emplace(spec.id, out.species.size());
        out.species.push_back(std::move(spec));
    }
    if (out.species.empty()) return fail("<model> declares no species");

    std::unordered_map<std::string, double> globals;
    for (const pugi::xml_node& p : children(child(model, "listOfParameters"), "parameter"))
        globals[p.attribute("id").as_string("")] = p.attribute("value").as_double(0.0);

    for (const pugi::xml_node& r : children(child(model, "listOfReactions"), "reaction")) {
        ReactionSpec rs;
        rs.id = r.attribute("id").as_string("");
        rs.reversible = r.attribute("reversible").as_bool(level == 2);
        auto side = [&](std::string_view list, std::vector<std::pair<std::string, double>>& into) {
            for (const pugi::xml_node& ref : children(child(r, list), "speciesReference")) {
                const std::string id = ref.attribute("species").as_string("");
                double stoich = ref.attribute("stoichiometry")
                                    ? ref.attribute("stoichiometry").as_double()
                                    : 1.0;
                if (child(ref, "stoichiometryMath")) {
                    if (error)
                        *error = "unsupported SBML construct <stoichiometryMath> on reaction '" +
                                 rs.id + "': a variable stoichiometry has no stoichiometric matrix";
                    stoich = std::nan("");
                }
                into.emplace_back(id, stoich);
            }
            return true;
        };
        side("listOfReactants", rs.reactants);
        side("listOfProducts", rs.products);
        for (const auto& [id, s] : rs.reactants)
            if (std::isnan(s)) return std::nullopt;
        for (const auto& [id, s] : rs.products)
            if (std::isnan(s)) return std::nullopt;
        for (const pugi::xml_node& m :
             children(child(r, "listOfModifiers"), "modifierSpeciesReference"))
            rs.modifiers.push_back(m.attribute("species").as_string(""));

        for (const auto& [id, s] : rs.reactants)
            if (!sym.speciesIndex.count(id))
                return fail("reaction '" + rs.id + "' references undeclared species '" + id + "'");
        for (const auto& [id, s] : rs.products)
            if (!sym.speciesIndex.count(id))
                return fail("reaction '" + rs.id + "' references undeclared species '" + id + "'");

        pugi::xml_node law = child(r, "kineticLaw");
        if (!law)
            return fail("reaction '" + rs.id +
                        "' has no <kineticLaw>: BioCAD integrates rates, and a reaction without "
                        "one cannot be simulated (only flux-balanced)");
        // Local parameters shadow globals, as SBML specifies.
        sym.parameters = globals;
        pugi::xml_node locals = child(law, "listOfLocalParameters");
        if (!locals) locals = child(law, "listOfParameters");
        for (const pugi::xml_node& p : children(locals, "localParameter"))
            sym.parameters[p.attribute("id").as_string("")] = p.attribute("value").as_double(0.0);
        for (const pugi::xml_node& p : children(locals, "parameter"))
            sym.parameters[p.attribute("id").as_string("")] = p.attribute("value").as_double(0.0);

        pugi::xml_node math = child(law, "math");
        if (!math) return fail("reaction '" + rs.id + "' has a <kineticLaw> with no <math>");
        pugi::xml_node expr;
        for (pugi::xml_node c : math.children())
            if (c.type() == pugi::node_element) { expr = c; break; }
        if (!expr) return fail("reaction '" + rs.id + "' has an empty <math>");

        // Refuse any operator outside the closed set the five laws use, naming it.
        static const char* kAllowed[] = {"times", "plus", "minus", "divide", "power"};
        std::vector<pugi::xml_node> stack{expr};
        while (!stack.empty()) {
            pugi::xml_node n = stack.back();
            stack.pop_back();
            const std::string_view op = applyOperator(n);
            if (!op.empty()) {
                bool ok = false;
                for (const char* a : kAllowed) ok = ok || op == a;
                if (!ok)
                    return fail("reaction '" + rs.id + "': kinetic law uses the MathML operator <" +
                                std::string(op) +
                                "/>, which is outside BioCAD's five rate laws (mass action, "
                                "reversible mass action, Michaelis-Menten, reversible "
                                "Michaelis-Menten, Hill)");
                for (const pugi::xml_node& a : applyArguments(n)) stack.push_back(a);
            }
        }

        const std::string substrate = rs.reactants.empty() ? std::string() : rs.reactants[0].first;
        const std::string productId = rs.products.empty() ? std::string() : rs.products[0].first;
        bool matched = false;

        // 1. reversible mass action: minus(product over reactants, product over products)
        if (applyOperator(expr) == "minus") {
            const auto args = applyArguments(expr);
            if (args.size() == 2) {
                const Product f = matchProduct(args[0], sym);
                const Product b = matchProduct(args[1], sym);
                if (f.ok && b.ok && powersMatch(f.powers, rs.reactants) &&
                    powersMatch(b.powers, rs.products)) {
                    rs.law = RateLaw::ReversibleMassAction;
                    rs.parameters = {f.constant, b.constant};
                    rs.parameterNames = {"kf", "kr"};
                    rs.reversible = true;
                    matched = true;
                }
            }
        }
        // 2. plain mass action
        if (!matched) {
            const Product p = matchProduct(expr, sym);
            if (p.ok && powersMatch(p.powers, rs.reactants)) {
                rs.law = RateLaw::MassAction;
                rs.parameters = {p.constant};
                rs.parameterNames = {"k"};
                matched = true;
            }
        }
        // 3. Hill before Michaelis-Menten: an exponent of 1 makes Hill degenerate to
        //    MM, and matching MM first would silently discard a fitted n.
        if (!matched && !substrate.empty()) {
            double vmax = 0, k = 0, nHill = 0;
            if (matchHill(expr, sym, substrate, vmax, k, nHill)) {
                rs.law = RateLaw::Hill;
                rs.parameters = {vmax, k, nHill};
                rs.parameterNames = {"Vmax", "K", "n"};
                matched = true;
            }
        }
        if (!matched && !substrate.empty()) {
            double vmax = 0, km = 0;
            if (matchMichaelisMenten(expr, sym, substrate, vmax, km)) {
                rs.law = RateLaw::MichaelisMenten;
                rs.parameters = {vmax, km};
                rs.parameterNames = {"Vmax", "Km"};
                matched = true;
            }
        }
        if (!matched && !substrate.empty() && !productId.empty()) {
            double vf = 0, ks = 0, vr = 0, kp = 0;
            if (matchReversibleMm(expr, sym, substrate, productId, vf, ks, vr, kp)) {
                rs.law = RateLaw::ReversibleMichaelisMenten;
                rs.parameters = {vf, ks, vr, kp};
                rs.parameterNames = {"Vf", "Kms", "Vr", "Kmp"};
                rs.reversible = true;
                matched = true;
            }
        }
        if (!matched)
            return fail("reaction '" + rs.id +
                        "': the kinetic law is well-formed MathML but is not one of BioCAD's five "
                        "rate laws - it may mix species that are not reactants, use an enzyme "
                        "modifier as a variable Vmax, or apply a compartment factor. It is "
                        "refused rather than approximated by the nearest supported form");
        out.reactions.push_back(std::move(rs));
    }
    if (out.reactions.empty()) return fail("<model> declares no reactions");
    return analyze(out);
}

std::string exportSbml(const NetworkSpec& network) {
    std::ostringstream o;
    o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version2/core\" level=\"3\" "
         "version=\"2\">\n";
    o << "  <model id=\"" << network.id << "\">\n";

    std::vector<std::string> compartments;
    for (const SpeciesSpec& s : network.species) {
        const std::string c = s.compartment.empty() ? std::string("default") : s.compartment;
        if (std::find(compartments.begin(), compartments.end(), c) == compartments.end())
            compartments.push_back(c);
    }
    o << "    <listOfCompartments>\n";
    for (const std::string& c : compartments)
        o << "      <compartment id=\"" << c << "\" size=\"1\" constant=\"true\"/>\n";
    o << "    </listOfCompartments>\n";

    o << "    <listOfSpecies>\n";
    for (const SpeciesSpec& s : network.species) {
        o << "      <species id=\"" << s.id << "\" name=\"" << s.name << "\" compartment=\""
          << (s.compartment.empty() ? "default" : s.compartment) << "\" initialConcentration=\""
          << fmt(s.initialConcentration) << "\" boundaryCondition=\""
          << (s.boundary ? "true" : "false") << "\" constant=\"false\" "
          << "hasOnlySubstanceUnits=\"false\"/>\n";
    }
    o << "    </listOfSpecies>\n";

    o << "    <listOfReactions>\n";
    for (const ReactionSpec& r : network.reactions) {
        o << "      <reaction id=\"" << r.id << "\" reversible=\""
          << (r.reversible ? "true" : "false") << "\">\n";
        auto side = [&](const char* tag, const std::vector<std::pair<std::string, double>>& v) {
            if (v.empty()) return;
            o << "        <listOf" << tag << ">\n";
            for (const auto& [id, s] : v)
                o << "          <speciesReference species=\"" << id << "\" stoichiometry=\""
                  << fmt(s) << "\" constant=\"true\"/>\n";
            o << "        </listOf" << tag << ">\n";
        };
        side("Reactants", r.reactants);
        side("Products", r.products);
        if (!r.modifiers.empty()) {
            o << "        <listOfModifiers>\n";
            for (const std::string& m : r.modifiers)
                o << "          <modifierSpeciesReference species=\"" << m << "\"/>\n";
            o << "        </listOfModifiers>\n";
        }
        // Parameter names are local to the reaction, so a network with two mass-action
        // steps does not need two globally unique names for "k".
        const char* const* names = rateLawParameterNames(r.law);
        o << "        <kineticLaw>\n"
          << "          <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
        auto product = [&](const char* indent, const std::string& constant,
                           const std::vector<std::pair<std::string, double>>& v) {
            o << indent << "<apply>\n" << indent << "  <times/>\n"
              << indent << "  <ci> " << constant << " </ci>\n";
            for (const auto& [id, s] : v) {
                if (std::abs(s - 1.0) < 1e-12) {
                    o << indent << "  <ci> " << id << " </ci>\n";
                } else {
                    o << indent << "  <apply>\n" << indent << "    <power/>\n"
                      << indent << "    <ci> " << id << " </ci>\n"
                      << indent << "    <cn> " << fmt(s) << " </cn>\n"
                      << indent << "  </apply>\n";
                }
            }
            o << indent << "</apply>\n";
        };
        switch (r.law) {
            case RateLaw::MassAction:
                product("            ", names[0], r.reactants);
                break;
            case RateLaw::ReversibleMassAction:
                o << "            <apply>\n              <minus/>\n";
                product("              ", names[0], r.reactants);
                product("              ", names[1], r.products);
                o << "            </apply>\n";
                break;
            case RateLaw::MichaelisMenten:
                o << "            <apply>\n              <divide/>\n";
                product("              ", names[0],
                        {{r.reactants.front().first, 1.0}});
                o << "              <apply>\n                <plus/>\n"
                  << "                <ci> " << names[1] << " </ci>\n"
                  << "                <ci> " << r.reactants.front().first << " </ci>\n"
                  << "              </apply>\n            </apply>\n";
                break;
            case RateLaw::Hill: {
                const std::string s = r.reactants.front().first;
                o << "            <apply>\n              <divide/>\n"
                  << "              <apply>\n                <times/>\n"
                  << "                <ci> " << names[0] << " </ci>\n"
                  << "                <apply>\n                  <power/>\n"
                  << "                  <ci> " << s << " </ci>\n"
                  << "                  <ci> " << names[2] << " </ci>\n"
                  << "                </apply>\n              </apply>\n"
                  << "              <apply>\n                <plus/>\n"
                  << "                <apply>\n                  <power/>\n"
                  << "                  <ci> " << names[1] << " </ci>\n"
                  << "                  <ci> " << names[2] << " </ci>\n"
                  << "                </apply>\n"
                  << "                <apply>\n                  <power/>\n"
                  << "                  <ci> " << s << " </ci>\n"
                  << "                  <ci> " << names[2] << " </ci>\n"
                  << "                </apply>\n              </apply>\n            </apply>\n";
                break;
            }
            case RateLaw::ReversibleMichaelisMenten: {
                const std::string s = r.reactants.front().first;
                const std::string p = r.products.front().first;
                o << "            <apply>\n              <divide/>\n"
                  << "              <apply>\n                <minus/>\n"
                  << "                <apply>\n                  <divide/>\n"
                  << "                  <apply>\n                    <times/>\n"
                  << "                    <ci> " << names[0] << " </ci>\n"
                  << "                    <ci> " << s << " </ci>\n"
                  << "                  </apply>\n"
                  << "                  <ci> " << names[1] << " </ci>\n"
                  << "                </apply>\n"
                  << "                <apply>\n                  <divide/>\n"
                  << "                  <apply>\n                    <times/>\n"
                  << "                    <ci> " << names[2] << " </ci>\n"
                  << "                    <ci> " << p << " </ci>\n"
                  << "                  </apply>\n"
                  << "                  <ci> " << names[3] << " </ci>\n"
                  << "                </apply>\n              </apply>\n"
                  << "              <apply>\n                <plus/>\n"
                  << "                <cn> 1 </cn>\n"
                  << "                <apply>\n                  <divide/>\n"
                  << "                  <ci> " << s << " </ci>\n"
                  << "                  <ci> " << names[1] << " </ci>\n"
                  << "                </apply>\n"
                  << "                <apply>\n                  <divide/>\n"
                  << "                  <ci> " << p << " </ci>\n"
                  << "                  <ci> " << names[3] << " </ci>\n"
                  << "                </apply>\n              </apply>\n            </apply>\n";
                break;
            }
        }
        o << "          </math>\n            <listOfLocalParameters>\n";
        for (std::size_t i = 0; i < r.parameters.size() && names[i]; ++i)
            o << "              <localParameter id=\"" << names[i] << "\" value=\""
              << fmt(r.parameters[i]) << "\"/>\n";
        o << "            </listOfLocalParameters>\n        </kineticLaw>\n";
        o << "      </reaction>\n";
    }
    o << "    </listOfReactions>\n  </model>\n</sbml>\n";
    return o.str();
}

}  // namespace biocad::sim
