// sim/Enrichment.h - pathway over-representation and network topology.
//
// THE BACKGROUND IS AN ARGUMENT, NOT A DEFAULT. A hypergeometric p-value is a
// function of four numbers, and one of them is the size of the universe the query was
// drawn from. Using "every gene in the database" when the experiment could only have
// detected 9000 transcripts inflates every p-value's significance, and it is the most
// common way an enrichment figure becomes wrong. enrich() therefore takes the
// background as a required argument, intersects the pathway sets with it, and records
// the resulting universe size in the report.
//
// BENJAMINI-HOCHBERG, NOT BONFERRONI, AND NEVER A RAW P-VALUE. A pathway database has
// hundreds to thousands of sets; a raw p of 0.01 over 1500 tests is an expectation of
// 15 false positives. The q-values here are the standard step-up BH values with
// monotonicity enforced, and the implementation is checked against the worked example
// in Benjamini & Hochberg's 1995 paper.
//
// GRAPH METRICS REQUIRE EVIDENCE PER EDGE. An edge with no stated evidence is DROPPED
// and named in GraphMetrics::warnings. Degree, betweenness and community membership
// are only as trustworthy as the edge set, and an interaction network assembled from
// unattributed edges produces confident numbers about nothing. No docking score, and
// nothing derived from one, may ever become an edge weight here.
#pragma once

#include <string>
#include <vector>

#include "data/Systems.h"

namespace biocad::sim {

// One gene-set pack in GMT form: id, description, then members, tab separated.
// `release` is stored on every hit, because a pathway result without the database
// version it came from cannot be reproduced.
struct GeneSetPack {
    std::string                                        release;
    std::vector<std::string>                           ids;
    std::vector<std::string>                           names;
    std::vector<std::vector<std::string>>              members;
    std::vector<std::string>                           warnings;
};

GeneSetPack parseGmt(const std::string& text, const std::string& release);

// Reads assets/packs/pathways/<file> through core::assetDir. `release` is taken from
// the pack's own leading comment lines (# release: ...) when present.
GeneSetPack loadGeneSetPack(const std::string& file);

// Upper-tail hypergeometric probability P(X >= k) for k successes drawn in `n` draws
// from a universe of `universe` containing `successes`. Computed with log-gammas, so
// a universe of 20000 does not overflow.
double hypergeometricUpperTail(int k, int n, int successes, int universe);

// Benjamini-Hochberg step-up q-values, in the input order of `pValues`.
std::vector<double> benjaminiHochberg(const std::vector<double>& pValues);

// The over-representation test. Query and background are de-duplicated; the query is
// intersected with the background (a query gene that is not in the background cannot
// be tested and is reported), and pathways are intersected with the background too.
EnrichmentReport enrich(const std::vector<std::string>& query,
                        const std::vector<std::string>& background, const GeneSetPack& pack,
                        int minimumPathwaySize = 3);

// Degree, connected components, Brandes betweenness and Louvain communities over an
// undirected view of the edge list.
GraphMetrics graph(const std::vector<NetworkEdge>& edges);

}  // namespace biocad::sim
