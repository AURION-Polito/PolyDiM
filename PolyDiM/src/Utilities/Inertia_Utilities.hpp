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
struct Inertia_Utilities final
{
    const Gedim::GeometryUtilities &geometryUtilities;

    struct Inertia_Data final
    {
        Eigen::Matrix3d Fmatrix;
        Eigen::Matrix3d FmatrixInv;
        Eigen::Vector3d translation;
        double absDetFmatrix;
        double signDetQ;
    };

    Inertia_Utilities(const Gedim::GeometryUtilities &geometryUtilities);

    void InertiaMapping2D(const Eigen::MatrixXd &vertices,
                          const Eigen::Vector3d &centroid,
                          const double &diameter,
                          const std::vector<Eigen::Matrix3d> &triangulation_vertices,
                          Inertia_Data &inertia_data) const;

    void InertiaMapping3D(const Eigen::MatrixXd &vertices,
                          const Eigen::Vector3d &centroid,
                          const double &diameter,
                          const std::vector<Eigen::MatrixXd> &tetrahedrons_vertices,
                          Inertia_Data &inertia_data) const;
};

} // namespace Utilities
} // namespace Polydim

#endif
