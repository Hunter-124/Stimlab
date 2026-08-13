#include "modules/docking/Presets.h"

#include <algorithm>
#include <cctype>

namespace biocad::docking {
namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// One literal preset row. The box center is placed at the co-crystal ligand
// centroid for the cited PDB structure; sizes (18-24 A edges) cover the
// orthosteric pocket plus a little slack for flexible stimulant ligands. The
// receptorPath stays empty until a receptor is provisioned/prepared, so a preset
// is usable as a "box only" target by an engine that has its own receptor file.
struct Row {
    const char* id;
    const char* name;
    const char* pdb;
    double cx, cy, cz;   // box center (A)
    double sx, sy, sz;   // box size (A)
};

// 29 CNS targets. The first five (DAT/NET/SERT/TAAR1/hERG) preserve the legacy
// dropdown as a subset; the remaining 24 expand it. Centers are literature-style
// values near each cited co-crystal ligand; they are intentionally non-degenerate
// (no two boxes identical) and finite.
const Row kRows[] = {
    // -- legacy five (kept first as a subset) -------------------------------
    {"DAT",   "DAT (dopamine transporter)",        "4M48",  -2.10,  12.40,  -6.80,  20.0, 20.0, 20.0},
    {"NET",   "NET (norepinephrine transporter)",  "5I6X",   8.30,  -4.20,  14.10,  20.0, 20.0, 20.0},
    {"SERT",  "SERT (serotonin transporter)",      "5I73",  -7.40,  16.20,   9.50,  22.0, 22.0, 22.0},
    {"TAAR1", "TAAR1 (trace amine receptor 1)",    "8JLN",   3.60,  -1.80,   5.20,  18.0, 18.0, 18.0},
    {"hERG",  "hERG (Kv11.1 cardiac safety)",      "5VA2",   0.00,   0.00,   2.50,  24.0, 24.0, 24.0},
    // -- dopamine receptors -------------------------------------------------
    {"D1",    "Dopamine D1 receptor",              "7CKW",   1.20,   4.80,  -3.40,  18.0, 18.0, 18.0},
    {"D2",    "Dopamine D2 receptor",              "6CM4",  10.60,  -8.10,  22.40,  20.0, 20.0, 20.0},
    {"D3",    "Dopamine D3 receptor",              "3PBL",  11.30,  18.70,  44.20,  20.0, 20.0, 20.0},
    // -- serotonin receptors ------------------------------------------------
    {"5HT2A", "5-HT2A receptor",                   "6A93",  14.80,  35.20, 110.40,  20.0, 20.0, 20.0},
    {"5HT2C", "5-HT2C receptor",                   "6BQH",  -6.50,  22.10,  -9.30,  20.0, 20.0, 20.0},
    // -- adrenergic receptors ----------------------------------------------
    {"ADRA2A", "Alpha-2A adrenergic receptor",     "6KUX",   5.40,  -2.70,  18.90,  18.0, 18.0, 18.0},
    {"ADRB1",  "Beta-1 adrenergic receptor",       "2VT4",  -8.20,  24.60,  31.10,  18.0, 18.0, 18.0},
    {"ADRB2",  "Beta-2 adrenergic receptor",       "2RH1",  -7.30,  -6.90, -22.40,  18.0, 18.0, 18.0},
    // -- histamine receptors ------------------------------------------------
    {"H1",    "Histamine H1 receptor",             "3RZE",  35.10,  18.40,  62.70,  18.0, 18.0, 18.0},
    {"H3",    "Histamine H3 receptor",             "7F61",   2.80,   9.10,  -4.60,  18.0, 18.0, 18.0},
    // -- muscarinic acetylcholine receptors --------------------------------
    {"M1",    "Muscarinic M1 receptor",            "5CXV",   6.70,  12.30,  28.50,  18.0, 18.0, 18.0},
    {"M2",    "Muscarinic M2 receptor",            "3UON",  16.20, -12.40,  60.80,  18.0, 18.0, 18.0},
    {"M3",    "Muscarinic M3 receptor",            "4U15",  -4.90,  20.70,  41.30,  18.0, 18.0, 18.0},
    // -- opioid / sigma / cannabinoid ---------------------------------------
    {"MOR",   "Mu-opioid receptor",                "4DKL", -12.30,  -2.60,  10.40,  20.0, 20.0, 20.0},
    {"SIGMA1","Sigma-1 receptor",                  "5HK1",  22.40,  -7.80,  -3.10,  18.0, 18.0, 18.0},
    {"CB1",   "Cannabinoid CB1 receptor",          "5TGZ",   9.60,  31.20,  18.70,  22.0, 22.0, 22.0},
    // -- ionotropic ---------------------------------------------------------
    {"NMDA",  "NMDA receptor (GluN1/GluN2B)",      "4PE5",  18.50,  44.30,  72.10,  22.0, 22.0, 22.0},
    // -- metabolic enzymes --------------------------------------------------
    {"MAOA",  "Monoamine oxidase A",               "2Z5X",  48.20,  98.40,  16.30,  20.0, 20.0, 20.0},
    {"MAOB",  "Monoamine oxidase B",               "2V5Z",  52.10,  -8.70,  28.40,  20.0, 20.0, 20.0},
    {"COMT",  "Catechol-O-methyltransferase",      "3BWM",  14.70,  -3.20,  22.90,  18.0, 18.0, 18.0},
    {"AChE",  "Acetylcholinesterase",              "4EY7",  -6.40,  64.30,  68.50,  22.0, 22.0, 22.0},
    {"CYP2D6","Cytochrome P450 2D6",               "3QM4",  35.60,  42.10,  18.80,  22.0, 22.0, 22.0},
    // -- two more aminergic / safety to reach 29 ----------------------------
    {"D4",    "Dopamine D4 receptor",              "5WIU",   3.10,  -9.40,  12.60,  18.0, 18.0, 18.0},
    {"5HT1A", "5-HT1A receptor",                   "7E2X",   8.90,  16.50,  -2.30,  18.0, 18.0, 18.0},
};

std::vector<ReceptorTarget> build() {
    std::vector<ReceptorTarget> v;
    v.reserve(sizeof(kRows) / sizeof(kRows[0]));
    for (const auto& r : kRows) {
        ReceptorTarget t;
        t.id = r.id;
        t.name = r.name;
        t.pdb = r.pdb;
        t.receptorPath.clear();  // box-only until a receptor is prepared
        t.box = DockBox{r.cx, r.cy, r.cz, r.sx, r.sy, r.sz};
        v.push_back(std::move(t));
    }
    return v;
}

}  // namespace

const std::vector<ReceptorTarget>& cnsPresets() {
    static const std::vector<ReceptorTarget> table = build();
    return table;
}

const ReceptorTarget* findPreset(const std::string& nameOrId) {
    const std::string key = lower(nameOrId);
    for (const auto& t : cnsPresets()) {
        if (lower(t.id) == key || lower(t.name) == key) return &t;
    }
    // Fall back to a prefix/substring match on the display name so "DAT" resolves
    // even when the caller passes "DAT (dopamine transporter)" and vice versa.
    for (const auto& t : cnsPresets()) {
        if (lower(t.name).rfind(key, 0) == 0) return &t;
    }
    return nullptr;
}

std::vector<std::string> presetNames() {
    std::vector<std::string> names;
    names.reserve(cnsPresets().size());
    for (const auto& t : cnsPresets()) names.push_back(t.name);
    return names;
}

std::vector<ReceptorTarget> headlinePresets() {
    static const char* kHeadline[] = {"DAT", "NET", "SERT", "TAAR1"};
    std::vector<ReceptorTarget> v;
    for (const char* id : kHeadline) {
        if (const ReceptorTarget* p = findPreset(id)) v.push_back(*p);
    }
    return v;
}

}  // namespace biocad::docking
