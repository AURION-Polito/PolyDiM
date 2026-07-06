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
/// @brief Compute the barycentric weights of the 1D Lagrange basis.
///
/// Returns the weights \f$w_i = \dfrac{1}{\prod_{j \neq i} (x_i - x_j)}\f$ associated
/// with the interpolation nodes \f$\{x_i\}\f$. These are computed once and reused by
/// Lagrange_1D_values() and Lagrange_1D_derivative_values() to evaluate the basis and
/// its derivatives efficiently.
///
/// @param interpolation_points_x Interpolation nodes \f$\{x_i\}\f$.
/// @return The vector of barycentric weights \f$w_i\f$ (empty if no nodes are given).
Eigen::VectorXd Lagrange_1D_coefficients(const Eigen::VectorXd &interpolation_points_x);

/// @brief Evaluate the 1D Lagrange basis functions at given points.
///
/// Computes \f$\ell_i(x) = w_i \prod_{j \neq i} (x - x_j)\f$ for every evaluation point
/// \f$x\f$ and every basis function \f$\ell_i\f$, where \f$w_i\f$ are the precomputed
/// barycentric weights. With a single interpolation node the basis reduces to the constant \f$1\f$.
///
/// @param interpolation_points_x   Interpolation nodes \f$\{x_i\}\f$.
/// @param lagrange_1D_coefficients Barycentric weights from Lagrange_1D_coefficients().
/// @param evaluation_points_x      Points at which to evaluate the basis.
/// @return A \f$N_{\text{eval}} \times N_{\text{nodes}}\f$ matrix whose \f$(p, i)\f$ entry
///         is \f$\ell_i\f$ evaluated at the \f$p\f$-th point.
Eigen::MatrixXd Lagrange_1D_values(const Eigen::VectorXd &interpolation_points_x,
                                   const Eigen::VectorXd &lagrange_1D_coefficients,
                                   const Eigen::VectorXd &evaluation_points_x);

/// @brief Evaluate the derivatives of the 1D Lagrange basis functions at given points.
///
/// Computes \f$\ell_i'(x) = w_i \sum_{j \neq i} \prod_{k \neq i,\,j} (x - x_k)\f$ (product
/// rule applied to Lagrange_1D_values()) for every evaluation point and every basis
/// function. With a single interpolation node the derivative is identically zero.
///
/// @param interpolation_points_x   Interpolation nodes \f$\{x_i\}\f$.
/// @param lagrange_1D_coefficients Barycentric weights from Lagrange_1D_coefficients().
/// @param evaluation_points_x      Points at which to evaluate the derivatives.
/// @return A \f$N_{\text{eval}} \times N_{\text{nodes}}\f$ matrix whose \f$(p, i)\f$ entry
///         is \f$\ell_i'\f$ evaluated at the \f$p\f$-th point.
Eigen::MatrixXd Lagrange_1D_derivative_values(const Eigen::VectorXd &interpolation_points_x,
                                              const Eigen::VectorXd &lagrange_1D_coefficients,
                                              const Eigen::VectorXd &evaluation_points_x);

} // namespace Lagrange
} // namespace Interpolation
} // namespace Polydim

#endif
