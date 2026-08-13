// bio/Rotamer.h - the backbone-dependent rotamer library, side-chain construction
// from internal coordinates, and Goldstein dead-end elimination.
//
// WHAT THIS SHIPS, AND WHAT IT DOES NOT. The Dunbrack 2010 smoothed
// backbone-dependent rotamer library is CC BY 4.0 and would be redistributable
// with attribution, but its only distribution channel is a download page issued
// after a manually approved licence application, so it could not be obtained for
// this tree. Inventing rotamer angles instead is not an option. The shipped pack
// (assets/packs/rotamers/rotamers-pdb-derived-2026.json) is therefore DERIVED: its
// means, circular standard deviations, probabilities, bond lengths and bond angles
// are measured from 219 X-ray entries at <= 1.4 A by scripts/build-rotamer-pack.py,
// and the pack says so in its own `notDunbrack` field. The binning convention
// follows Shapovalov & Dunbrack 2011 (Structure 19:844-858); the numbers are not
// theirs and are not interchangeable with theirs. docs/variants.md states exactly
// what is and is not included.
//
// THERE IS NO ENERGY HERE. The selection objective below is a unit-free ranking
// built from two things only: heavy-atom clash counts and the library probability
// of a rotamer. It is never reported as a number - RotamerRebuild has no energy
// field, deliberately - and it must not be turned into one.
#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "bio/Structure.h"

namespace biocad::bio {

// ---------------------------------------------------------------------------
// The library
// ---------------------------------------------------------------------------

struct RotamerEntry {
    std::string         name;          // "mt", "p", "mttt", ...
    int                 count = 0;     // observations behind it
    double              probability = 0;
    std::vector<double> chi;           // degrees, chi1..chiN
    std::vector<double> chiSd;         // circular standard deviation, degrees
};

// One side-chain atom's internal coordinates: place it at `bondLength` from
// parents[2], making `bondAngle` with parents[1]-parents[2], at torsion
// parents[0]-parents[1]-parents[2]-atom. `chi` is 0 when that torsion is the
// measured mean (a rigid part: a ring atom, a branch, an improper) and k > 0 when
// the torsion IS chi_k and therefore varies with the rotamer.
struct BuildAtom {
    std::string atom;
    std::string parents[3];
    double      bondLength = 0;
    double      bondAngle = 0;
    double      torsion = 0;
    double      torsionSd = 0;
    int         chi = 0;
    int         observations = 0;
};

struct PhiPsiBin {
    int phi = 0;      // lower edge, degrees
    int psi = 0;
    int count = 0;
    std::vector<RotamerEntry> rotamers;
};

struct ResidueRotamers {
    std::string name;
    int chiCount = 0;
    int nonRotamericChi = 0;   // index of the continuous terminal chi, 0 if none
    int symmetricChi = 0;
    std::vector<BuildAtom>    build;
    std::vector<RotamerEntry> backboneIndependent;
    int                       backboneIndependentCount = 0;
    std::vector<PhiPsiBin>    bins;
};

class RotamerLibrary {
public:
    [[nodiscard]] bool empty() const { return residues_.empty(); }
    [[nodiscard]] const std::string& name() const { return name_; }
    [[nodiscard]] const std::string& method() const { return method_; }
    [[nodiscard]] const std::vector<std::string>& attribution() const { return attribution_; }
    [[nodiscard]] int binDegrees() const { return binDegrees_; }
    [[nodiscard]] int datasetEntryCount() const { return datasetEntries_; }
    [[nodiscard]] int residuesMeasured() const { return residuesMeasured_; }

    [[nodiscard]] const ResidueRotamers* residue(const std::string& threeLetter) const;

    // Rotamers for `threeLetter` at this backbone. When the phi/psi bin has too
    // few observations to have been emitted, the backbone-independent aggregate is
    // returned and `usedBackboneIndependent` is set - which the caller reports as
    // an assumption rather than silently absorbing.
    [[nodiscard]] std::vector<RotamerEntry> rotamersAt(const std::string& threeLetter,
                                                       double phi, double psi,
                                                       bool& usedBackboneIndependent) const;

    friend RotamerLibrary parseRotamerLibrary(const nlohmann::json& j);

private:
    std::string name_, method_;
    std::vector<std::string> attribution_;
    int binDegrees_ = 60;
    int datasetEntries_ = 0;
    int residuesMeasured_ = 0;
    std::vector<ResidueRotamers> residues_;
};

// Throws core::Error on a wrong schemaVersion or a structurally invalid document.
RotamerLibrary parseRotamerLibrary(const nlohmann::json& j);
RotamerLibrary loadRotamerLibrary(const std::filesystem::path& file);

// The pack shipped beside the executable. Returns nullptr when it is missing or
// malformed: a missing library is a refused rebuild naming the file, never a
// side chain built from guessed angles.
const RotamerLibrary* builtinRotamerLibrary();

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

struct Point3d {
    double x = 0, y = 0, z = 0;
};

// NeRF placement: the point at `bond` from c, angle b-c-d = `angleDeg`, torsion
// a-b-c-d = `torsionDeg` (IUPAC sign).
[[nodiscard]] Point3d placeAtom(const Point3d& a, const Point3d& b, const Point3d& c,
                                double bond, double angleDeg, double torsionDeg);

[[nodiscard]] double dihedralDegrees(const Point3d& a, const Point3d& b, const Point3d& c,
                                     const Point3d& d);

// phi (C(-1)-N-CA-C) and psi (N-CA-C-N(+1)) for residue `i` of a chain. Returns
// false when a terminal residue leaves one of them undefined - which is a real
// case, not an error, and the caller falls back to the backbone-independent set.
[[nodiscard]] bool backboneTorsions(const Chain& chain, std::size_t index, double& phi,
                                    double& psi);

// Builds the side chain of `threeLetter` on the N/CA/C backbone of `backbone`,
// with the given chi angles in degrees. Returns the built heavy atoms (CB
// onwards). Empty when the residue has no build tree in the library or the
// backbone is incomplete.
[[nodiscard]] std::vector<Atom> buildSideChain(const ResidueRotamers& tmpl,
                                               const Residue& backbone,
                                               const std::vector<double>& chi);

// Two heavy atoms clash when they are closer than `overlapFactor` times the sum of
// their van der Waals radii (bio::vdwRadius). 0.75 is the factor used throughout
// this module: at 0.75 a C-C pair clashes below 2.55 A, which is well inside any
// real non-bonded contact and outside every bonded 1-3 distance the builder emits.
inline constexpr double kClashOverlapFactor = 0.75;

[[nodiscard]] int countClashes(const std::vector<Atom>& a, const std::vector<Atom>& b,
                               double overlapFactor = kClashOverlapFactor);

// ---------------------------------------------------------------------------
// Dead-end elimination
// ---------------------------------------------------------------------------

// A packing problem in the ONLY currency this module has: `self[i][r]` is
// position i / rotamer r's own cost (clashes against everything held fixed, plus
// the library-probability term), and `pair[i][j][r][s]` is the clash count between
// two rotamers at two flexible positions. Both are unit-free by construction and
// neither is reported.
struct DeeProblem {
    std::vector<std::vector<double>> self;
    // pair[i][j][r][s], filled for every ordered (i, j) with i != j and symmetric:
    // pair[i][j][r][s] == pair[j][i][s][r].
    std::vector<std::vector<std::vector<std::vector<double>>>> pair;

    [[nodiscard]] std::size_t positions() const { return self.size(); }
    [[nodiscard]] double total(const std::vector<int>& choice) const;
};

struct DeeResult {
    // alive[i][r] is 0 once rotamer r at position i has been eliminated.
    std::vector<std::vector<char>> alive;
    int  eliminated = 0;
    int  passes = 0;
    // The chosen rotamer per position after the final repack over the survivors.
    std::vector<int> chosen;
    double totalScore = 0;
    bool exhaustive = false;   // false when the survivor space forced iterative repack
};

// Goldstein 1994 (Biophys J 66:1335-1340): rotamer r at position i is eliminated
// when some competitor t at the same position is better no matter what the other
// positions do, i.e.
//
//   self[i][r] - self[i][t] + sum_{j != i} min_s ( pair[i][j][r][s] - pair[i][j][t][s] ) > 0
//
// The sum is over live rotamers s only, and elimination repeats until a pass
// removes nothing. This is exact: an eliminated rotamer cannot be in the global
// optimum. The final repack is exhaustive when the surviving space is at most
// `exhaustiveLimit` combinations, and otherwise iterative single-position descent
// from the best-self start, which DeeResult::exhaustive reports.
[[nodiscard]] DeeResult goldsteinDee(const DeeProblem& problem,
                                     std::size_t exhaustiveLimit = 200000);

}  // namespace biocad::bio
