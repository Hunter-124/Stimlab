// bio/Superpose.h - Kabsch optimal rigid-body superposition and RMSD.
//
// Pure geometry over paired point sets. It deliberately knows nothing about
// bio::Structure or chem::Conformer: the caller selects which atoms correspond
// (e.g. CA atoms of aligned residues) and passes coordinates. There is no
// overload taking a chem::Conformer, because a distance-geometry embedding of a
// small molecule is not a protein and must never reach a structure-comparison
// number.
#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace biocad::bio {

using Point3 = std::array<double, 3>;

// Rotation is row-major 3x3. Applying the fit to a mobile point p:
//   q = rotation * (p - mobileCentroid) + referenceCentroid
// `translation` is the equivalent single vector for the unshifted form
//   q = rotation * p + translation
struct Superposition {
    std::array<double, 9> rotation{1, 0, 0, 0, 1, 0, 0, 0, 1};
    Point3 translation{0, 0, 0};
    Point3 mobileCentroid{0, 0, 0};
    Point3 referenceCentroid{0, 0, 0};
    double rmsd = 0.0;          // Angstrom, over the paired points
    std::size_t pairs = 0;
    // True when the raw SVD product was a reflection and the third column of V
    // had to be negated. Kabsch without this correction returns a mirror image
    // with a deceptively low RMSD, which is a wrong answer, not a warning.
    bool reflectionCorrected = false;
};

// Throws core::Error on size mismatch or fewer than 3 points (fewer than three
// points does not determine a rotation).
Superposition kabsch(const std::vector<Point3>& mobile, const std::vector<Point3>& reference);

// Applies a computed fit to a point.
Point3 applySuperposition(const Superposition& s, const Point3& p);

// RMSD of two point sets in their given frames, with NO superposition. Useful
// for scoring an externally produced alignment.
double rmsdInPlace(const std::vector<Point3>& a, const std::vector<Point3>& b);

}  // namespace biocad::bio
