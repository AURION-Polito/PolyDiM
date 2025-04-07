#ifndef __VEM_Inertia_Utilities_HPP
#define __VEM_Inertia_Utilities_HPP

#include "GeometryUtilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
struct VEM_Inertia_Utilities final
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

    VEM_Inertia_Utilities(const Gedim::GeometryUtilities &geometryUtilities) : geometryUtilities(geometryUtilities)
    {
    }

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
} // namespace VEM
} // namespace Polydim

#endif
