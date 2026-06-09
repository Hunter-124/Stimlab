#include "modules/docking/ReceptorPrep.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <unordered_set>
#include <utility>

#include "chem/Molecule.h"
#include "core/AppPaths.h"
#include "modules/docking/EngineLocator.h"
#include "modules/docking/PdbqtWriter.h"

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace stimlab::docking {
namespace {

namespace fs = std::filesystem;

// Ring atom names that are aromatic carbons/nitrogens in the standard aromatic
// residues. Used so PHE/TYR/TRP/HIS ring carbons get the AutoDock 'A' (aromatic C)
// type and the imidazole/indole nitrogens get 'NA'.
bool isAromaticRingAtom(const std::string& res, const std::string& atom) {
    static const std::set<std::pair<std::string, std::string>> ar = {
        {"PHE", "CG"}, {"PHE", "CD1"}, {"PHE", "CD2"}, {"PHE", "CE1"}, {"PHE", "CE2"}, {"PHE", "CZ"},
        {"TYR", "CG"}, {"TYR", "CD1"}, {"TYR", "CD2"}, {"TYR", "CE1"}, {"TYR", "CE2"}, {"TYR", "CZ"},
        {"TRP", "CG"}, {"TRP", "CD1"}, {"TRP", "CD2"}, {"TRP", "NE1"}, {"TRP", "CE2"}, {"TRP", "CE3"},
        {"TRP", "CZ2"}, {"TRP", "CZ3"}, {"TRP", "CH2"},
        {"HIS", "CG"}, {"HIS", "ND1"}, {"HIS", "CD2"}, {"HIS", "CE1"}, {"HIS", "NE2"},
        {"HID", "CG"}, {"HID", "ND1"}, {"HID", "CD2"}, {"HID", "CE1"}, {"HID", "NE2"},
        {"HIE", "CG"}, {"HIE", "ND1"}, {"HIE", "CD2"}, {"HIE", "CE1"}, {"HIE", "NE2"},
        {"HIP", "CG"}, {"HIP", "ND1"}, {"HIP", "CD2"}, {"HIP", "CE1"}, {"HIP", "NE2"},
    };
    return ar.count({res, atom}) > 0;
}

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

bool isWater(const std::string& res) {
    return res == "HOH" || res == "WAT" || res == "DOD" || res == "H2O" || res == "TIP" ||
           res == "SOL" || res == "TIP3" || res == "TIP4";
}

// Infer atomic number from the PDB element field (cols 77-78) if present/valid,
// else from the atom name. Receptor atoms we keep are protein heavy atoms
// (C/N/O/S/P) so single-letter inference from the atom name is reliable.
int inferElement(const std::string& elemField, const std::string& atomName) {
    const std::string e = upper(trim(elemField));
    if (!e.empty()) {
        // Normalise case to the element table's convention before lookup.
        std::string sym = e;
        if (sym.size() == 2) sym[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(sym[1])));
        if (const auto* info = chem::findElementBySymbol(sym)) return info->z;
    }
    // Fall back to the atom name: first alphabetic character is the element for the
    // protein subset (e.g. "CA"->C alpha carbon, "OG1"->O, "ND1"->N, "SD"->S).
    for (char c : atomName) {
        const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        switch (u) {
            case 'C': return 6;
            case 'N': return 7;
            case 'O': return 8;
            case 'S': return 16;
            case 'P': return 15;
            case 'H': return 1;
            default: break;
        }
        if (std::isalpha(static_cast<unsigned char>(c))) break;  // unknown leading element
    }
    return 6;  // carbon-like default; never 0/NaN
}

// Approximate, element-based partial charge. Vina's scoring is charge-independent
// so this column is informational; values are kept small and finite for AutoDock4-
// style consumers. (A converged protein PEOE solve is out of scope without RDKit.)
double approxReceptorCharge(int z, bool aromatic) {
    switch (z) {
        case 8:  return -0.30;             // O (OA)
        case 7:  return aromatic ? -0.20 : -0.30;  // N / NA
        case 16: return -0.10;             // S / SA
        case 6:  return aromatic ? 0.00 : 0.08;    // C / A
        case 15: return 0.30;              // P
        case 1:  return 0.20;              // polar H (only if a tool added them)
        default: return 0.00;
    }
}

double parseCol(const std::string& line, size_t start, size_t len) {
    if (line.size() < start) return std::nan("");
    const std::string s = trim(line.substr(start, std::min(len, line.size() - start)));
    if (s.empty()) return std::nan("");
    try {
        return std::stod(s);
    } catch (...) {
        return std::nan("");
    }
}

std::string colStr(const std::string& line, size_t start, size_t len) {
    if (line.size() <= start) return {};
    return trim(line.substr(start, std::min(len, line.size() - start)));
}

bool startsWith(const std::string& s, const char* p) { return s.rfind(p, 0) == 0; }

#if defined(_WIN32)
// Best-effort PowerShell file download (mirrors EngineLocator::ensureVina). Returns
// true iff a non-trivial file now exists at `dest`. Never throws; quoted paths.
bool downloadFile(const std::string& url, const fs::path& dest, std::string& note) {
    std::string ps;
    ps += "$ErrorActionPreference='Stop';";
    ps += "try{";
    ps += "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;";
    ps += "Invoke-WebRequest -UseBasicParsing -Uri '" + url + "' -OutFile '" + dest.string() + "';";
    ps += "exit 0}catch{exit 2}";
    const std::string cmd = "powershell -NoProfile -NonInteractive -Command \"" + ps + "\"";
    const int code = std::system(cmd.c_str());
    std::error_code ec;
    if (code == 0 && fs::exists(dest, ec) && fs::file_size(dest, ec) > 256) {
        return true;
    }
    note = "download failed (exit " + std::to_string(code) + ")";
    if (fs::exists(dest, ec) && fs::file_size(dest, ec) <= 256) fs::remove(dest, ec);  // error page
    return false;
}
#endif

std::string readFile(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// The structure-frame docking-box remark the docking module reads back to override
// the preset center. Empty when no box could be derived.
std::string boxRemark(bool has, double x, double y, double z) {
    if (!has) return {};
    char b[96];
    std::snprintf(b, sizeof(b), "REMARK STIMLAB_BOX %.3f %.3f %.3f\n", x, y, z);
    return b;
}

}  // namespace

fs::path receptorsDir() {
    return AppPaths::instance().runtime() / "receptors";
}

ReceptorPdbqt pdbToRigidReceptor(const std::string& pdbText,
                                 const std::vector<std::string>& keepHetero) {
    ReceptorPdbqt out;
    std::unordered_set<std::string> keep;
    for (const auto& h : keepHetero) keep.insert(upper(h));

    std::istringstream in(pdbText);
    std::string line;
    std::string body;
    body.reserve(pdbText.size());
    int serial = 0;
    int modelSeen = 0;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (startsWith(line, "MODEL")) {
            if (++modelSeen > 1) break;  // NMR ensemble: keep only the first model
            continue;
        }
        if (startsWith(line, "ENDMDL")) {
            if (modelSeen >= 1) break;
            continue;
        }
        if (startsWith(line, "END")) break;

        const bool isAtom = startsWith(line, "ATOM");
        const bool isHet = startsWith(line, "HETATM");
        if (!isAtom && !isHet) continue;

        // Fixed PDB columns (0-based): name 12-15, altLoc 16, resName 17-19,
        // chain 21, resSeq 22-25, x 30-37, y 38-45, z 46-53, element 76-77.
        const std::string atomName = colStr(line, 12, 4);
        const char altLoc = line.size() > 16 ? line[16] : ' ';
        const std::string resName = upper(colStr(line, 17, 3));
        char chain = line.size() > 21 ? line[21] : 'A';
        if (chain == ' ') chain = 'A';
        const std::string resSeq = colStr(line, 22, 4);

        // Keep a single conformation: blank / 'A' / '1' alternate location only.
        if (!(altLoc == ' ' || altLoc == 'A' || altLoc == '1')) continue;

        // Drop ALL HETATM (waters, cryo/buffer additives, AND the co-crystal ligand
        // that defined the box) unless explicitly retained as a cofactor. An empty
        // keep set therefore yields a clean apo-like pocket. An explicit keep wins.
        if (isHet && !keep.count(resName)) {
            ++out.droppedHetero;
            continue;
        }

        const double x = parseCol(line, 30, 8);
        const double y = parseCol(line, 38, 8);
        const double z = parseCol(line, 46, 8);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;

        const std::string elemField = line.size() > 76 ? line.substr(76, 2) : std::string();
        const int zNum = inferElement(elemField, atomName);
        if (zNum == 1) continue;  // drop hydrogens (built-in path is heavy-atom)

        const bool aromatic = isAromaticRingAtom(resName, atomName);
        const std::string adType = autodockAtomType(zNum, aromatic, /*polarH=*/false);
        const double charge = approxReceptorCharge(zNum, aromatic);

        int resNum = 0;
        try { resNum = resSeq.empty() ? 0 : std::stoi(resSeq); } catch (...) { resNum = 0; }

        ++serial;
        char buf[128];
        // Fixed-width PDBQT ATOM line: standard PDB columns 1-66 preserved, then the
        // PDBQT trailing partial charge (%10.3f) + right-justified 2-char atom type.
        std::snprintf(buf, sizeof(buf),
                      "ATOM  %5d %-4.4s %-3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10.3f %-2.2s\n",
                      serial, atomName.c_str(), resName.c_str(), chain, resNum, x, y, z,
                      1.00, 0.00, charge, adType.c_str());
        body += buf;

        ++out.atomCount;
        if (isHet) ++out.keptHetero;
    }

    std::string text;
    text.reserve(body.size() + 256);
    text += "REMARK  StimLab rigid receptor (built-in heavy-atom prep; AutoDock4 types)\n";
    text += "REMARK  binding-affinity / target-engagement prediction only; not a synthesis artifact\n";
    text += "REMARK  Vina scoring is partial-charge-independent; charges are approximate\n";
    text += body;
    out.text = std::move(text);

    out.note = "kept " + std::to_string(out.atomCount) + " receptor atoms";
    if (out.keptHetero) out.note += " (incl. " + std::to_string(out.keptHetero) + " cofactor)";
    out.note += ", dropped " + std::to_string(out.droppedHetero) + " water/hetero atoms";

    // Derive a docking box center in THIS structure's coordinate frame so the box
    // overlaps the receptor regardless of the preset's literature-style center.
    out.hasBox = receptorBoxFromPdb(pdbText, out.cx, out.cy, out.cz, out.boxSource);
    if (out.hasBox) out.note += "; box from " + out.boxSource;
    return out;
}

bool receptorBoxFromPdb(const std::string& pdbText, double& cx, double& cy, double& cz,
                        std::string& source) {
    struct Acc { int n = 0; double sx = 0, sy = 0, sz = 0; };
    std::map<std::string, Acc> hetGroups;  // candidate co-crystal ligands by residue
    Acc protein;
    int modelSeen = 0;

    std::istringstream in(pdbText);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (startsWith(line, "MODEL")) { if (++modelSeen > 1) break; continue; }
        if (startsWith(line, "ENDMDL")) { if (modelSeen >= 1) break; continue; }
        if (startsWith(line, "END")) break;

        const bool isAtom = startsWith(line, "ATOM");
        const bool isHet = startsWith(line, "HETATM");
        if (!isAtom && !isHet) continue;

        const char altLoc = line.size() > 16 ? line[16] : ' ';
        if (!(altLoc == ' ' || altLoc == 'A' || altLoc == '1')) continue;
        const std::string resName = upper(colStr(line, 17, 3));
        const double x = parseCol(line, 30, 8), y = parseCol(line, 38, 8), z = parseCol(line, 46, 8);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;

        if (isAtom) {
            ++protein.n; protein.sx += x; protein.sy += y; protein.sz += z;
        } else if (!isWater(resName)) {
            const std::string chain = line.size() > 21 ? std::string(1, line[21]) : "A";
            const std::string key = resName + "/" + chain + "/" + colStr(line, 22, 4);
            Acc& a = hetGroups[key];
            ++a.n; a.sx += x; a.sy += y; a.sz += z;
        }
    }

    // Prefer the largest non-water HETATM residue with >= 8 atoms (the co-crystal
    // ligand) - its centroid is the canonical binding-site center.
    const Acc* best = nullptr;
    std::string bestKey;
    for (const auto& [key, a] : hetGroups) {
        if (a.n >= 8 && (!best || a.n > best->n)) { best = &a; bestKey = key; }
    }
    if (best) {
        cx = best->sx / best->n; cy = best->sy / best->n; cz = best->sz / best->n;
        source = "co-crystal ligand " + bestKey.substr(0, bestKey.find('/'));
        return true;
    }
    if (protein.n > 0) {
        cx = protein.sx / protein.n; cy = protein.sy / protein.n; cz = protein.sz / protein.n;
        source = "receptor centroid";
        return true;
    }
    return false;
}

ReceptorPrepResult locatePreparedReceptor(const std::string& targetId) {
    ReceptorPrepResult r;
    std::error_code ec;
    const fs::path pdbqt = receptorsDir() / (targetId + ".pdbqt");
    if (fs::exists(pdbqt, ec) && fs::is_regular_file(pdbqt, ec) && fs::file_size(pdbqt, ec) > 0) {
        r.ready = true;
        r.path = pdbqt.string();
        r.note = "prepared receptor cached at " + r.path;
        // Recover the structure-frame box center from the STIMLAB_BOX remark.
        std::ifstream f(pdbqt);
        std::string l;
        while (std::getline(f, l)) {
            if (l.rfind("REMARK STIMLAB_BOX", 0) == 0) {
                std::istringstream is(l.substr(std::string("REMARK STIMLAB_BOX").size()));
                if (is >> r.cx >> r.cy >> r.cz) r.hasBox = true;
                break;
            }
            if (l.rfind("ATOM", 0) == 0 || l.rfind("HETATM", 0) == 0) break;  // remarks are first
        }
    } else {
        r.note = "no prepared receptor for " + targetId + " (provision to fetch + prepare)";
    }
    return r;
}

ReceptorPrepResult ensureReceptor(const ReceptorTarget& target, bool allowDownload) {
    // 1) Cached prepared receptor wins (no network).
    ReceptorPrepResult r = locatePreparedReceptor(target.id);
    if (r.ready) return r;

    if (!allowDownload) return r;  // locate-only: stay in fallback with the note above

    std::error_code ec;
    const fs::path dir = receptorsDir();
    fs::create_directories(dir, ec);
    if (ec) {
        r.note = "could not create receptors dir: " + ec.message();
        return r;
    }
    const fs::path pdbqtPath = dir / (target.id + ".pdbqt");
    const fs::path pdbPath = dir / (target.pdb + ".pdb");

#if defined(_WIN32)
    // 2) Fetch the source PDB from RCSB if not already cached.
    if (!(fs::exists(pdbPath, ec) && fs::file_size(pdbPath, ec) > 256)) {
        if (target.pdb.empty()) {
            r.note = "preset " + target.id + " has no PDB id to fetch.";
            return r;
        }
        const std::string url = "https://files.rcsb.org/download/" + target.pdb + ".pdb";
        std::string dlNote;
        if (!downloadFile(url, pdbPath, dlNote)) {
            r.note = "PDB " + target.pdb + " fetch failed (" + dlNote +
                     "); staying in descriptor-estimate fallback.";
            return r;
        }
    }
#else
    if (!(fs::exists(pdbPath, ec) && fs::file_size(pdbPath, ec) > 256)) {
        r.note = "PDB fetch only implemented on Windows.";
        return r;
    }
#endif

    // 3) Read the structure once + derive the box center in its own coordinate frame.
    const std::string pdbText = readFile(pdbPath);
    if (pdbText.empty()) {
        r.note = "fetched PDB " + target.pdb + " was empty/unreadable.";
        return r;
    }
    double bx = 0, by = 0, bz = 0;
    std::string boxSrc;
    const bool hasBox = receptorBoxFromPdb(pdbText, bx, by, bz, boxSrc);
    const std::string boxLine = boxRemark(hasBox, bx, by, bz);

#if defined(_WIN32)
    // 4a) Preferred path: obabel.exe is provisioned -> higher-quality prep with H.
    if (auto obabel = locateEngine(Engine::Obabel)) {
        const std::wstring cmd = L"\"" + obabel->wstring() + L"\" \"" + pdbPath.wstring() +
                                 L"\" -xr -O \"" + pdbqtPath.wstring() + L"\"";
        STARTUPINFOW si{};
        si.cb = sizeof(si);
        PROCESS_INFORMATION pi{};
        std::wstring buf = cmd;
        if (CreateProcessW(nullptr, buf.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr,
                           dir.c_str(), &si, &pi)) {
            WaitForSingleObject(pi.hProcess, 120000);
            CloseHandle(pi.hProcess);
            CloseHandle(pi.hThread);
        }
        if (fs::exists(pdbqtPath, ec) && fs::file_size(pdbqtPath, ec) > 0) {
            // Prepend our box remark so the dock path can recover the center.
            const std::string obOut = readFile(pdbqtPath);
            std::ofstream o(pdbqtPath, std::ios::binary | std::ios::trunc);
            o << boxLine << obOut;
            o.close();
            r.ready = true;
            r.path = pdbqtPath.string();
            r.usedObabel = true;
            r.cx = bx; r.cy = by; r.cz = bz; r.hasBox = hasBox;
            r.note = "prepared " + target.id + " from " + target.pdb + " via obabel (rigid receptor)" +
                     (hasBox ? "; box from " + boxSrc : "") + ".";
            return r;
        }
        // obabel present but produced nothing usable -> fall through to built-in.
    }
#endif

    // 4b) Built-in path: parse the cached PDB -> rigid receptor PDBQT.
    ReceptorPdbqt rec = pdbToRigidReceptor(pdbText);
    if (rec.empty()) {
        r.note = "built-in prep parsed no receptor atoms from " + target.pdb + ".";
        return r;
    }
    {
        std::ofstream o(pdbqtPath, std::ios::binary);
        o << boxLine << rec.text;
    }
    if (!(fs::exists(pdbqtPath, ec) && fs::file_size(pdbqtPath, ec) > 0)) {
        r.note = "failed to write prepared receptor for " + target.id + ".";
        return r;
    }
    r.ready = true;
    r.path = pdbqtPath.string();
    r.atomCount = rec.atomCount;
    r.cx = rec.cx; r.cy = rec.cy; r.cz = rec.cz; r.hasBox = rec.hasBox;
    r.note = "prepared " + target.id + " from " + target.pdb + " (built-in heavy-atom writer; " +
             rec.note + ").";
    return r;
}

}  // namespace stimlab::docking
