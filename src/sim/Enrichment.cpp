#include "sim/Enrichment.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "core/Assets.h"

namespace biocad::sim {
namespace {

double logChoose(int n, int k) {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

std::vector<std::string> unique(const std::vector<std::string>& v) {
    std::vector<std::string> out = v;
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

}  // namespace

double hypergeometricUpperTail(int k, int n, int successes, int universe) {
    if (universe <= 0 || n <= 0 || successes <= 0) return 1.0;
    if (k <= 0) return 1.0;
    const int top = std::min(n, successes);
    if (k > top) return 0.0;
    const double denom = logChoose(universe, n);
    double sum = 0;
    for (int i = k; i <= top; ++i) {
        const double term = logChoose(successes, i) + logChoose(universe - successes, n - i) - denom;
        if (std::isfinite(term)) sum += std::exp(term);
    }
    return std::min(1.0, std::max(0.0, sum));
}

std::vector<double> benjaminiHochberg(const std::vector<double>& pValues) {
    const std::size_t m = pValues.size();
    std::vector<double> q(m, 1.0);
    if (m == 0) return q;
    std::vector<std::size_t> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return pValues[a] < pValues[b]; });
    // Step-up from the largest p, carrying the running minimum: that is what makes
    // the q-values monotone in p, which a naive p*m/rank is not.
    double running = 1.0;
    for (std::size_t rank = m; rank >= 1; --rank) {
        const std::size_t idx = order[rank - 1];
        const double scaled = pValues[idx] * static_cast<double>(m) / static_cast<double>(rank);
        running = std::min(running, scaled);
        q[idx] = std::min(1.0, running);
    }
    return q;
}

GeneSetPack parseGmt(const std::string& text, const std::string& release) {
    GeneSetPack pack;
    pack.release = release;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            const std::string tag = "# release:";
            if (line.rfind(tag, 0) == 0) {
                std::string v = line.substr(tag.size());
                const std::size_t a = v.find_first_not_of(" \t");
                if (a != std::string::npos) pack.release = v.substr(a);
            }
            continue;
        }
        std::vector<std::string> fields;
        std::string field;
        std::istringstream ls(line);
        while (std::getline(ls, field, '\t')) fields.push_back(field);
        if (fields.size() < 3) {
            pack.warnings.push_back("skipped a GMT line with fewer than three fields: '" +
                                    line.substr(0, 40) + "'");
            continue;
        }
        // Broad's GMT convention is id-then-name; Reactome ships name-then-id. The
        // stable id is unmistakable (R-HSA-164843), so detect it rather than making
        // every caller declare which dialect its file is in.
        bool reactomeOrder = fields[1].rfind("R-", 0) == 0 &&
                             fields[1].find('-', 2) != std::string::npos &&
                             fields[1].size() > 6;
        pack.ids.push_back(reactomeOrder ? fields[1] : fields[0]);
        pack.names.push_back(reactomeOrder ? fields[0] : fields[1]);
        pack.members.push_back(
            unique(std::vector<std::string>(fields.begin() + 2, fields.end())));
    }
    return pack;
}

GeneSetPack loadGeneSetPack(const std::string& file) {
    const auto dir = core::assetDir("packs/pathways");
    std::ifstream in(dir / file);
    if (!in) {
        GeneSetPack pack;
        pack.warnings.push_back("gene-set pack '" + file + "' not found under " + dir.string());
        return pack;
    }
    std::ostringstream o;
    o << in.rdbuf();
    return parseGmt(o.str(), "unstated (the pack carried no '# release:' line)");
}

EnrichmentReport enrich(const std::vector<std::string>& query,
                        const std::vector<std::string>& background, const GeneSetPack& pack,
                        int minimumPathwaySize) {
    EnrichmentReport out;
    out.method = "hypergeometric (upper tail) + Benjamini-Hochberg";
    out.databaseRelease = pack.release;
    out.background = unique(background);
    const std::vector<std::string> q = unique(query);

    if (out.background.empty()) {
        out.warnings.push_back(
            "no background given: the hypergeometric answer is a function of the universe the "
            "query was drawn from, so there is no result to report");
        return out;
    }
    std::unordered_set<std::string> bg(out.background.begin(), out.background.end());
    std::vector<std::string> tested, untestable;
    for (const std::string& g : q) {
        if (bg.count(g)) tested.push_back(g);
        else untestable.push_back(g);
    }
    out.querySet = tested;
    if (!untestable.empty()) {
        std::string list;
        for (std::size_t i = 0; i < untestable.size() && i < 8; ++i)
            list += (i ? ", " : "") + untestable[i];
        out.warnings.push_back(std::to_string(untestable.size()) +
                               " query identifier(s) are not in the background and cannot be "
                               "tested (" + list +
                               (untestable.size() > 8 ? ", ..." : "") + ")");
    }
    if (tested.empty()) {
        out.warnings.push_back("no query identifier survived intersection with the background");
        return out;
    }
    std::unordered_set<std::string> querySet(tested.begin(), tested.end());
    const int universe = static_cast<int>(out.background.size());
    const int drawn = static_cast<int>(tested.size());

    std::vector<EnrichmentHit> hits;
    std::vector<double> ps;
    int skippedSmall = 0;
    for (std::size_t i = 0; i < pack.ids.size(); ++i) {
        // The pathway is intersected with the background too: a pathway member the
        // experiment could not have detected is not a member for this test.
        int size = 0, overlap = 0;
        for (const std::string& g : pack.members[i]) {
            if (!bg.count(g)) continue;
            ++size;
            if (querySet.count(g)) ++overlap;
        }
        if (size < minimumPathwaySize) { ++skippedSmall; continue; }
        if (overlap == 0) continue;
        EnrichmentHit h;
        h.pathwayId = pack.ids[i];
        h.pathwayName = pack.names[i];
        h.inSetAndPathway = overlap;
        h.pathwaySize = size;
        h.pValue = hypergeometricUpperTail(overlap, drawn, size, universe);
        const double expected = static_cast<double>(drawn) * size / universe;
        h.foldEnrichment = expected > 0 ? overlap / expected : 0.0;
        h.source = pack.release;
        hits.push_back(std::move(h));
        ps.push_back(hits.back().pValue);
    }
    const std::vector<double> qs = benjaminiHochberg(ps);
    for (std::size_t i = 0; i < hits.size(); ++i) hits[i].qValue = qs[i];
    std::sort(hits.begin(), hits.end(), [](const EnrichmentHit& a, const EnrichmentHit& b) {
        if (a.qValue != b.qValue) return a.qValue < b.qValue;
        return a.pValue < b.pValue;
    });
    out.hits = std::move(hits);
    // The number of tests IS part of the result: the q-values cannot be recomputed
    // without it.
    out.warnings.push_back("tested " + std::to_string(ps.size()) + " pathway(s) against a " +
                           std::to_string(universe) + "-identifier background with " +
                           std::to_string(drawn) + " query identifiers; " +
                           std::to_string(skippedSmall) +
                           " pathway(s) fell below the minimum size of " +
                           std::to_string(minimumPathwaySize) + " after intersection");
    return out;
}

GraphMetrics graph(const std::vector<NetworkEdge>& edges) {
    GraphMetrics out;
    std::vector<NetworkEdge> kept;
    for (const NetworkEdge& e : edges) {
        if (e.evidence.empty() || e.provenance == Provenance::NotComputed) {
            out.warnings.push_back("dropped the edge " + e.source + "-" + e.target +
                                   ": an edge with no stated evidence cannot support a degree, a "
                                   "betweenness or a community assignment");
            continue;
        }
        if (e.source == e.target) {
            out.warnings.push_back("dropped the self-loop on " + e.source +
                                   ": it changes no shortest path and inflates the degree");
            continue;
        }
        kept.push_back(e);
    }
    out.edges = kept;

    std::unordered_map<std::string, std::size_t> index;
    auto id = [&](const std::string& name) {
        auto it = index.find(name);
        if (it != index.end()) return it->second;
        index.emplace(name, out.nodes.size());
        out.nodes.push_back(name);
        return out.nodes.size() - 1;
    };
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    std::vector<double> weights;
    for (const NetworkEdge& e : kept) {
        const std::size_t a = id(e.source), b = id(e.target);
        pairs.emplace_back(a, b);
        weights.push_back(e.weight != 0 ? e.weight : 1.0);
    }
    const std::size_t n = out.nodes.size();
    if (n == 0) return out;

    std::vector<std::vector<std::size_t>> adjacency(n);
    std::vector<std::vector<double>> adjacencyWeight(n);
    out.degree.assign(n, 0);
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        const auto [a, b] = pairs[k];
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
        adjacencyWeight[a].push_back(weights[k]);
        adjacencyWeight[b].push_back(weights[k]);
        ++out.degree[a];
        ++out.degree[b];
    }

    // Connected components by breadth-first sweep.
    std::vector<int> component(n, -1);
    int components = 0;
    for (std::size_t s = 0; s < n; ++s) {
        if (component[s] != -1) continue;
        std::queue<std::size_t> queue;
        queue.push(s);
        component[s] = components;
        while (!queue.empty()) {
            const std::size_t v = queue.front();
            queue.pop();
            for (std::size_t w : adjacency[v])
                if (component[w] == -1) { component[w] = components; queue.push(w); }
        }
        ++components;
    }
    out.componentCount = components;

    // Brandes (2001) betweenness on the unweighted graph. Halved at the end because
    // the accumulation counts every unordered pair twice on an undirected graph.
    out.betweenness.assign(n, 0.0);
    for (std::size_t s = 0; s < n; ++s) {
        std::vector<std::vector<std::size_t>> predecessors(n);
        std::vector<double> sigma(n, 0.0);
        std::vector<long long> distance(n, -1);
        std::vector<std::size_t> stack;
        std::queue<std::size_t> queue;
        sigma[s] = 1.0;
        distance[s] = 0;
        queue.push(s);
        while (!queue.empty()) {
            const std::size_t v = queue.front();
            queue.pop();
            stack.push_back(v);
            for (std::size_t w : adjacency[v]) {
                if (distance[w] < 0) {
                    distance[w] = distance[v] + 1;
                    queue.push(w);
                }
                if (distance[w] == distance[v] + 1) {
                    sigma[w] += sigma[v];
                    predecessors[w].push_back(v);
                }
            }
        }
        std::vector<double> delta(n, 0.0);
        for (std::size_t i = stack.size(); i-- > 0;) {
            const std::size_t w = stack[i];
            for (std::size_t v : predecessors[w])
                delta[v] += sigma[v] / sigma[w] * (1.0 + delta[w]);
            if (w != s) out.betweenness[w] += delta[w];
        }
    }
    for (double& b : out.betweenness) b *= 0.5;

    // Louvain: local moving to maximise modularity, then aggregate, repeated until no
    // node moves. Weighted, because a confidence-weighted interaction network is the
    // normal case.
    double totalWeight = 0;
    for (double w : weights) totalWeight += w;
    out.community.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) out.community[i] = static_cast<int>(i);
    if (totalWeight > 0) {
        std::vector<double> strength(n, 0.0);
        for (std::size_t v = 0; v < n; ++v)
            for (double w : adjacencyWeight[v]) strength[v] += w;
        const double m2 = 2.0 * totalWeight;
        std::vector<double> communityStrength(n, 0.0);
        for (std::size_t v = 0; v < n; ++v) communityStrength[v] = strength[v];
        bool moved = true;
        for (int pass = 0; pass < 100 && moved; ++pass) {
            moved = false;
            for (std::size_t v = 0; v < n; ++v) {
                const int from = out.community[v];
                communityStrength[from] -= strength[v];
                std::unordered_map<int, double> linkTo;
                for (std::size_t k = 0; k < adjacency[v].size(); ++k)
                    linkTo[out.community[adjacency[v][k]]] += adjacencyWeight[v][k];
                int best = from;
                double bestGain = linkTo.count(from) ? linkTo[from] -
                                      communityStrength[from] * strength[v] / m2
                                                     : -communityStrength[from] * strength[v] / m2;
                for (const auto& [c, link] : linkTo) {
                    const double gain = link - communityStrength[c] * strength[v] / m2;
                    if (gain > bestGain + 1e-12) { bestGain = gain; best = c; }
                }
                out.community[v] = best;
                communityStrength[best] += strength[v];
                if (best != from) moved = true;
            }
        }
        // Renumber the surviving communities to 0..k-1 so the labels mean something.
        std::unordered_map<int, int> relabel;
        for (int& c : out.community) {
            auto it = relabel.find(c);
            if (it == relabel.end()) it = relabel.emplace(c, static_cast<int>(relabel.size())).first;
            c = it->second;
        }
        // Modularity of the final partition, so the caller can judge it rather than
        // trusting the label.
        double q = 0;
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const auto [a, b] = pairs[k];
            if (out.community[a] == out.community[b]) q += 2.0 * weights[k];
        }
        std::unordered_map<int, double> degreeSum;
        for (std::size_t v = 0; v < n; ++v) degreeSum[out.community[v]] += strength[v];
        double expected = 0;
        for (const auto& [c, d] : degreeSum) expected += d * d;
        out.modularity = q / m2 - expected / (m2 * m2);
    }
    return out;
}

}  // namespace biocad::sim
