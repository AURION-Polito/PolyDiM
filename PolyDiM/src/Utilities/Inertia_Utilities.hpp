// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __Inertia_Utilities_HPP
#define __Inertia_Utilities_HPP

#include "GeometryUtilities.hpp"

namespace Polydim
{
namespace Utilities
{

/// @brief Affine map data aligning a polytope with its principal axes of inertia.
///
/// Encodes the affine change of coordinates
/// \f$\boldsymbol{x} = F\,\hat{\boldsymbol{x}} + \boldsymbol{t}\f$ (physical
/// \f$\leftrightarrow\f$ reference/inertia frame), where the reference frame is
/// isotropic and aligned with the element's principal axes of inertia. The mapping
/// is used to improve the conditioning of the monomial basis in inertia-based local spaces.
struct Inertia_Data final
{
    Eigen::Matrix3d Fmatrix;
    Eigen::Matrix3d FmatrixInv;
    Eigen::Vector3d translation;
    double absDetFmatrix;
    double signDetQ;
};

struct Inertia_Utilities final
{
    /// @brief Compute the inertia mapping of a polygon (2D).
    ///
    /// Constructs the affine map that sends the polygon onto an isotropic reference
    /// configuration aligned with its principal axes of inertia. The element is first
    /// rescaled by its diameter and translated to the centroid; the polygon mass
    /// (second-moment) matrix is then assembled from the triangulation and
    /// diagonalized (self-adjoint eigensolver), and its eigenpairs define a whitening
    /// transform that normalizes and aligns the principal directions; a final diameter
    /// rescaling is applied. The resulting map, its inverse, translation, Jacobian
    /// magnitude and orientation sign are stored in @p inertia_data.
    ///
    /// @param geometryUtilities      GeDiM helper providing polygon mass and diameter.
    /// @param vertices               Polygon vertices (one per column).
    /// @param centroid               Polygon centroid.
    /// @param diameter               Polygon diameter \f$h_E\f$ (initial rescaling factor).
    /// @param triangulation_vertices Sub-triangulation of the polygon (per-triangle vertices),
    ///                               used to compute the mass matrix.
    /// @param[out] inertia_data      Resulting affine-map data (see Inertia_Data).
    /// @throws std::runtime_error if the inertia (mass) matrix is singular.
    void InertiaMapping2D(const Gedim::GeometryUtilities &geometryUtilities,
                          const Eigen::MatrixXd &vertices,
                          const Eigen::Vector3d &centroid,
                          const double &diameter,
                          const std::vector<Eigen::Matrix3d> &triangulation_vertices,
                          Polydim::Utilities::Inertia_Data &inertia_data) const;

    /// @brief Compute the inertia mapping of a polyhedron (3D).
    ///
    /// Three-dimensional counterpart of InertiaMapping2D(): builds the affine map that
    /// aligns the polyhedron with its principal axes of inertia and normalizes it to an
    /// isotropic reference configuration, starting from a tetrahedral decomposition used
    /// to assemble the mass (second-moment) matrix. The map, its inverse, translation,
    /// Jacobian magnitude and orientation sign are stored in @p inertia_data.
    ///
    /// @param geometryUtilities     GeDiM helper providing polyhedron mass and diameter.
    /// @param vertices              Polyhedron vertices (one per column).
    /// @param centroid              Polyhedron centroid.
    /// @param diameter              Polyhedron diameter \f$h_E\f$ (initial rescaling factor).
    /// @param tetrahedrons_vertices Tetrahedral decomposition of the polyhedron (per-tetrahedron vertices),
    ///                              used to compute the mass matrix.
    /// @param[out] inertia_data     Resulting affine-map data (see Inertia_Data).
    /// @throws std::runtime_error if the inertia (mass) matrix is singular.
    void InertiaMapping3D(const Gedim::GeometryUtilities &geometryUtilities,
                          const Eigen::MatrixXd &vertices,
                          const Eigen::Vector3d &centroid,
                          const double &diameter,
                          const std::vector<Eigen::MatrixXd> &tetrahedrons_vertices,
                          Polydim::Utilities::Inertia_Data &inertia_data) const;
};

} // namespace Utilities
} // namespace Polydim

#endif
