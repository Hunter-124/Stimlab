// tests/test_enrichment.cpp - the hypergeometric test, Benjamini-Hochberg against the
// worked example in the original paper, the required background, and the graph metrics.
//
// Brandes betweenness is checked against a five-node path graph whose values can be
// counted by hand (0, 3, 4, 3, 0), because a betweenness implementation that is merely
// plausible is indistinguishable from a correct one on a real network.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "sim/Enrichment.h"

using namespace biocad;
using namespace biocad::sim;

TEST_CASE("Benjamini-Hochberg reproduces the 1995 paper's worked example",
          "[enrichment]") {
    // The 15 p-values of Needleman et al. as tabulated in Benjamini & Hochberg (1995),
    // J R Statist Soc B 57:289-300. At alpha = 0.05 the procedure rejects the four
    // smallest: i = 4 is the largest index with p_i <= (i/15)*0.05 (0.0095 <= 0.01333),
    // even though p_5 = 0.0201 > 0.01667.
    const std::vector<double> p = {0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298,
                                  0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590,
                                  1.0000};
    const std::vector<double> q = benjaminiHochberg(p);
    REQUIRE(q.size() == p.size());
    int rejected = 0;
    for (double v : q)
        if (v <= 0.05) ++rejected;
    REQUIRE(rejected == 4);
    // The step-up values, hand-computed as p_i * 15 / i with the running minimum.
    REQUIRE(std::abs(q[0] - 0.0015) < 1e-12);
    REQUIRE(std::abs(q[1] - 0.0030) < 1e-12);
    REQUIRE(std::abs(q[2] - 0.0095) < 1e-12);
    REQUIRE(std::abs(q[3] - 0.035625) < 1e-12);
    REQUIRE(std::abs(q[4] - 0.0603) < 1e-12);
    REQUIRE(std::abs(q[14] - 1.0) < 1e-12);
    // Monotone in p, which a naive p*m/rank is not: the running minimum from p_8 pulls
    // rank 7's value down from 0.0639 to 0.063857.
    for (std::size_t i = 1; i < q.size(); ++i) REQUIRE(q[i] >= q[i - 1] - 1e-15);
    REQUIRE(q[5] == q[6]);
    REQUIRE(benjaminiHochberg({}).empty());
}

TEST_CASE("the hypergeometric upper tail is exact on a countable case", "[enrichment]") {
    // Universe 10, 5 successes, 3 draws:
    //   P(X >= 2) = (C(5,2)C(5,1) + C(5,3)) / C(10,3) = (50 + 10)/120 = 0.5
    REQUIRE(std::abs(hypergeometricUpperTail(2, 3, 5, 10) - 0.5) < 1e-12);
    REQUIRE(std::abs(hypergeometricUpperTail(1, 3, 5, 10) - (1.0 - 10.0 / 120.0)) < 1e-12);
    REQUIRE(std::abs(hypergeometricUpperTail(3, 3, 5, 10) - 10.0 / 120.0) < 1e-12);
    // More successes than draws is impossible, not merely unlikely.
    REQUIRE(hypergeometricUpperTail(4, 3, 5, 10) == 0.0);
    // A 20000-gene universe must not overflow the binomials.
    const double big = hypergeometricUpperTail(5, 50, 100, 20000);
    REQUIRE(big > 0.0);
    REQUIRE(big < 1.0);
}

TEST_CASE("enrichment intersects with the background and refuses to run without one",
          "[enrichment]") {
    const std::string gmt =
        "# release: test fixture\n"
        "P1\tGlycolysis\tHK1\tPFKL\tALDOA\tGAPDH\tPKM\tENO1\n"
        "P2\tTCA cycle\tCS\tACO2\tIDH2\tOGDH\tSDHA\tFH\n"
        "P3\tTiny set\tAAA\tBBB\n";
    const GeneSetPack pack = parseGmt(gmt, "test fixture");
    REQUIRE(pack.ids.size() == 3);
    REQUIRE(pack.release == "test fixture");

    std::vector<std::string> background;
    for (int i = 0; i < 200; ++i) background.push_back("BG" + std::to_string(i));
    for (const auto& set : pack.members)
        for (const std::string& g : set) background.push_back(g);
    const std::vector<std::string> query = {"HK1",  "PFKL", "ALDOA",
                                            "GAPDH", "BG7",  "BG9", "NOT_IN_BACKGROUND"};

    const EnrichmentReport report = enrich(query, background, pack);
    REQUIRE(report.databaseRelease == "test fixture");
    REQUIRE_FALSE(report.hits.empty());
    REQUIRE(report.hits.front().pathwayId == "P1");
    // The 2-member set falls below the minimum size after intersection and is not
    // tested at all, so it cannot contribute a test to the BH correction.
    for (const EnrichmentHit& h : report.hits) REQUIRE(h.pathwayId != "P3");
    // The p-value must be the hypergeometric of the ACTUAL four numbers: universe 214,
    // 6 drawn (the untestable identifier is excluded), pathway size 6, overlap 4.
    REQUIRE(report.background.size() == 214);
    REQUIRE(report.querySet.size() == 6);
    REQUIRE(std::abs(report.hits.front().pValue - hypergeometricUpperTail(4, 6, 6, 214)) < 1e-15);
    // The untestable identifier is reported, not silently dropped.
    bool named = false;
    for (const std::string& w : report.warnings)
        if (w.find("NOT_IN_BACKGROUND") != std::string::npos) named = true;
    REQUIRE(named);

    const EnrichmentReport none = enrich(query, {}, pack);
    REQUIRE(none.hits.empty());
    REQUIRE_FALSE(none.warnings.empty());
}

TEST_CASE("Brandes betweenness matches a hand-counted five-node graph", "[enrichment]") {
    auto edge = [](const char* a, const char* b) {
        NetworkEdge e;
        e.source = a;
        e.target = b;
        e.weight = 1.0;
        e.evidence = "test fixture: an explicit adjacency";
        e.provenance = Provenance::Measured;
        return e;
    };
    // Path graph n1-n2-n3-n4-n5. Counting each unordered pair once:
    //   n2 lies on (1,3) (1,4) (1,5)             -> 3
    //   n3 lies on (1,4) (1,5) (2,4) (2,5)       -> 4
    //   n4 lies on (1,5) (2,5) (3,5)             -> 3
    //   n1 and n5 lie on nothing                 -> 0
    const std::vector<NetworkEdge> path = {edge("n1", "n2"), edge("n2", "n3"), edge("n3", "n4"),
                                           edge("n4", "n5")};
    const GraphMetrics g = graph(path);
    REQUIRE(g.nodes.size() == 5);
    REQUIRE(g.componentCount == 1);
    const double expected[5] = {0, 3, 4, 3, 0};
    for (std::size_t i = 0; i < 5; ++i) REQUIRE(std::abs(g.betweenness[i] - expected[i]) < 1e-12);
    REQUIRE(g.degree[0] == 1);
    REQUIRE(g.degree[2] == 2);

    // Two disjoint triangles: two components, two Louvain communities, modularity
    // exactly 0.5, and every betweenness zero because no vertex of a triangle is a cut
    // point.
    const std::vector<NetworkEdge> triangles = {edge("a", "b"), edge("b", "c"), edge("c", "a"),
                                                edge("x", "y"), edge("y", "z"), edge("z", "x")};
    const GraphMetrics t = graph(triangles);
    REQUIRE(t.componentCount == 2);
    for (double b : t.betweenness) REQUIRE(std::abs(b) < 1e-12);
    REQUIRE(t.community[0] == t.community[1]);
    REQUIRE(t.community[1] == t.community[2]);
    REQUIRE(t.community[3] == t.community[4]);
    REQUIRE(t.community[4] == t.community[5]);
    REQUIRE(t.community[0] != t.community[3]);
    REQUIRE(std::abs(t.modularity - 0.5) < 1e-9);
}

TEST_CASE("an edge with no evidence is dropped and named", "[enrichment]") {
    NetworkEdge withEvidence;
    withEvidence.source = "a";
    withEvidence.target = "b";
    withEvidence.weight = 1.0;
    withEvidence.evidence = "test fixture";
    withEvidence.provenance = Provenance::Measured;

    NetworkEdge bare = withEvidence;
    bare.source = "q";
    bare.target = "r";
    bare.evidence.clear();

    NetworkEdge unprovenanced = withEvidence;
    unprovenanced.source = "s";
    unprovenanced.target = "t";
    unprovenanced.provenance = Provenance::NotComputed;

    NetworkEdge loop = withEvidence;
    loop.source = "a";
    loop.target = "a";

    const GraphMetrics g = graph({withEvidence, bare, unprovenanced, loop});
    // Only a and b survive: an unattributed edge cannot support a degree.
    REQUIRE(g.nodes.size() == 2);
    REQUIRE(g.edges.size() == 1);
    REQUIRE(g.warnings.size() == 3);
    REQUIRE(g.degree[0] == 1);
}

TEST_CASE("the shipped Reactome pack loads and finds glycolysis", "[enrichment]") {
    const GeneSetPack pack = loadGeneSetPack("reactome-human.gmt");
    REQUIRE(pack.warnings.empty());
    // Release 97 ships 2868 human pathways; the assertion is deliberately loose on the
    // count and strict on the things a wrong parse would break.
    REQUIRE(pack.ids.size() > 2000);
    REQUIRE(pack.release.find("Reactome release") != std::string::npos);
    // Reactome's GMT is name-then-id, the reverse of the Broad convention, and the
    // loader must detect that rather than storing the name as the id.
    REQUIRE(pack.ids.front().rfind("R-HSA-", 0) == 0);
    REQUIRE_FALSE(pack.names.front().empty());

    std::vector<std::string> background;
    for (const auto& set : pack.members)
        for (const std::string& g : set) background.push_back(g);
    const std::vector<std::string> glycolysis = {"HK1",  "HK2",   "GPI",   "PFKL",
                                                "ALDOA", "TPI1", "GAPDH", "PGK1",
                                                "PGAM1", "ENO1", "PKM",   "LDHA"};
    const EnrichmentReport r = enrich(glycolysis, background, pack);
    REQUIRE_FALSE(r.hits.empty());
    REQUIRE(r.hits.front().pathwayName.find("Glycolysis") != std::string::npos);
    REQUIRE(r.hits.front().qValue < 1e-6);
    REQUIRE(r.hits.front().foldEnrichment > 10.0);
    REQUIRE(r.method == "hypergeometric (upper tail) + Benjamini-Hochberg");
}
