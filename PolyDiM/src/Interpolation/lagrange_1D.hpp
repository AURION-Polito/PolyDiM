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

#ifndef __Interpolation_Lagrange_1D_HPP
#define __Interpolation_Lagrange_1D_HPP

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
namespace Interpolation
{
namespace Lagrange
{
Eigen::VectorXd Lagrange_1D_coefficients(const Eigen::VectorXd &interpolation_points_x);
Eigen::MatrixXd Lagrange_1D_values(const Eigen::VectorXd &interpolation_points_x,
                                   const Eigen::VectorXd &lagrange_1D_coefficients,
                                   const Eigen::VectorXd &evaluation_points_x);
Eigen::MatrixXd Lagrange_1D_derivative_values(const Eigen::VectorXd &interpolation_points_x,
                                              const Eigen::VectorXd &lagrange_1D_coefficients,
                                              const Eigen::VectorXd &evaluation_points_x);

} // namespace Lagrange
} // namespace Interpolation
} // namespace Polydim

#endif
