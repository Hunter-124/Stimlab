#include "bio/Superpose.h"

#include <cmath>

#include <Eigen/Dense>

#include "core/Error.h"

namespace biocad::bio {
namespace {

Eigen::Vector3d centroidOf(const std::vector<Point3>& p) {
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (const auto& q : p) c += Eigen::Vector3d(q[0], q[1], q[2]);
    return c / static_cast<double>(p.size());
}

}  // namespace

Superposition kabsch(const std::vector<Point3>& mobile, const std::vector<Point3>& reference) {
    if (mobile.size() != reference.size()) {
        throw Error::invalidArgument("kabsch: point sets differ in size");
    }
    if (mobile.size() < 3) {
        throw Error::invalidArgument("kabsch: need at least 3 paired points to fit a rotation");
    }

    const std::size_t n = mobile.size();
    const Eigen::Vector3d cm = centroidOf(mobile);
    const Eigen::Vector3d cr = centroidOf(reference);

    // Centred coordinates, 3 x n so the covariance is a plain 3x3 product.
    Eigen::Matrix3Xd P(3, n), Q(3, n);
    for (std::size_t i = 0; i < n; ++i) {
        P.col(static_cast<Eigen::Index>(i)) =
            Eigen::Vector3d(mobile[i][0], mobile[i][1], mobile[i][2]) - cm;
        Q.col(static_cast<Eigen::Index>(i)) =
            Eigen::Vector3d(reference[i][0], reference[i][1], reference[i][2]) - cr;
    }

    // H = P * Q^T; R = V * U^T with H = U S V^T maximises trace(R * H), i.e.
    // minimises the RMSD (Kabsch 1976).
    const Eigen::Matrix3d H = P * Q.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    // The reflection correction. Without it the SVD solution may be an improper
    // rotation (det = -1), which mirrors the mobile set and can report a much
    // lower RMSD than any real rigid motion can achieve - a silently wrong
    // answer. Negating the third column of V (the direction of least variance)
    // is the minimal repair.
    bool corrected = false;
    if ((V * U.transpose()).determinant() < 0.0) {
        V.col(2) = -V.col(2);
        corrected = true;
    }
    const Eigen::Matrix3d R = V * U.transpose();

    Superposition out;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) out.rotation[static_cast<std::size_t>(r * 3 + c)] = R(r, c);
    }
    const Eigen::Vector3d t = cr - R * cm;
    out.translation = {t.x(), t.y(), t.z()};
    out.mobileCentroid = {cm.x(), cm.y(), cm.z()};
    out.referenceCentroid = {cr.x(), cr.y(), cr.z()};
    out.pairs = n;
    out.reflectionCorrected = corrected;

    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d d =
            R * P.col(static_cast<Eigen::Index>(i)) - Q.col(static_cast<Eigen::Index>(i));
        sum += d.squaredNorm();
    }
    out.rmsd = std::sqrt(sum / static_cast<double>(n));
    return out;
}

Point3 applySuperposition(const Superposition& s, const Point3& p) {
    const auto& r = s.rotation;
    return {r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + s.translation[0],
            r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + s.translation[1],
            r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + s.translation[2]};
}

double rmsdInPlace(const std::vector<Point3>& a, const std::vector<Point3>& b) {
    if (a.size() != b.size()) throw Error::invalidArgument("rmsdInPlace: size mismatch");
    if (a.empty()) throw Error::invalidArgument("rmsdInPlace: empty point sets");
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double dx = a[i][0] - b[i][0];
        const double dy = a[i][1] - b[i][1];
        const double dz = a[i][2] - b[i][2];
        sum += dx * dx + dy * dy + dz * dz;
    }
    return std::sqrt(sum / static_cast<double>(a.size()));
}

}  // namespace biocad::bio
