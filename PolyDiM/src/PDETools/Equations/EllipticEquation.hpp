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

#ifndef __PDETOOLS_EQUATION_EllipticEquation_HPP
#define __PDETOOLS_EQUATION_EllipticEquation_HPP

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
namespace PDETools
{
namespace Equations
{
struct EllipticEquation final
{
    Eigen::MatrixXd ComputeCellDiffusionMatrix(const Eigen::VectorXd &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &trial_basis_functions_derivative_values,
                                               const std::vector<Eigen::MatrixXd> &test_basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        Eigen::MatrixXd cell_matrix = test_basis_functions_derivative_values.at(0).transpose() *
                                      quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                                      trial_basis_functions_derivative_values.at(0);

        for (unsigned int d = 1; d < test_basis_functions_derivative_values.size(); ++d)
        {
            cell_matrix.noalias() += test_basis_functions_derivative_values.at(d).transpose() *
                                     quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                                     trial_basis_functions_derivative_values.at(d);
        }

        return cell_matrix;
    }

    Eigen::MatrixXd ComputeCellDiffusionMatrix(const Eigen::VectorXd &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        Eigen::MatrixXd cell_matrix = basis_functions_derivative_values.at(0).transpose() *
                                      quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                                      basis_functions_derivative_values.at(0);

        for (unsigned int d = 1; d < basis_functions_derivative_values.size(); ++d)
        {
            cell_matrix.noalias() += basis_functions_derivative_values.at(d).transpose() *
                                     quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                                     basis_functions_derivative_values.at(d);
        }

        return cell_matrix;
    }

    Eigen::MatrixXd ComputeCellDiffusionMatrix(const std::array<Eigen::VectorXd, 9> &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        const unsigned int numBasisFunctions = basis_functions_derivative_values[0].cols();
        const unsigned int dimension = basis_functions_derivative_values.size();

        Eigen::MatrixXd cell_matrix = Eigen::MatrixXd::Zero(numBasisFunctions, numBasisFunctions);
        for (unsigned int d1 = 0; d1 < dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < dimension; d2++)
            {
                cell_matrix.noalias() += basis_functions_derivative_values.at(d1).transpose() *
                                         quadrature_weights.cwiseProduct(diffusion_term_values.at(d1 + 3 * d2)).asDiagonal() *
                                         basis_functions_derivative_values.at(d2);
            }
        }
        return cell_matrix;
    }

    inline Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd &reaction_term_values,
                                                     const Eigen::MatrixXd &basis_functions_values,
                                                     const Eigen::VectorXd &quadrature_weights) const
    {
        return basis_functions_values.transpose() * quadrature_weights.cwiseProduct(reaction_term_values).asDiagonal() * basis_functions_values;
    }

    inline Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd &reaction_term_values,
                                                     const Eigen::MatrixXd &trial_basis_functions_values,
                                                     const Eigen::MatrixXd &test_basis_functions_values,
                                                     const Eigen::VectorXd &quadrature_weights) const
    {
        return test_basis_functions_values.transpose() *
            quadrature_weights.cwiseProduct(reaction_term_values).asDiagonal() *
            trial_basis_functions_values;
    }

    Eigen::MatrixXd ComputeCellAdvectionMatrix(const std::array<Eigen::VectorXd, 3> &advection_term_values,
                                               const Eigen::MatrixXd &test_basis_functions_values,
                                               const std::vector<Eigen::MatrixXd> &trial_basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        Eigen::MatrixXd cell_matrix = test_basis_functions_values.transpose() *
                                      quadrature_weights.cwiseProduct(advection_term_values[0]).asDiagonal() *
                                      trial_basis_functions_derivative_values[0];

        for (unsigned int d = 1; d < trial_basis_functions_derivative_values.size(); ++d)
        {
            cell_matrix.noalias() += test_basis_functions_values.transpose() *
                                     quadrature_weights.cwiseProduct(advection_term_values.at(d)).asDiagonal() *
                                     trial_basis_functions_derivative_values.at(d);
        }

        return cell_matrix;
    }

    inline Eigen::VectorXd ComputeCellForcingTerm(const Eigen::VectorXd &forcing_term_values,
                                                  const Eigen::MatrixXd &test_basis_functions_values,
                                                  const Eigen::VectorXd &quadrature_weights) const
    {
        return test_basis_functions_values.transpose() * quadrature_weights.asDiagonal() * forcing_term_values;
    }

    inline Eigen::VectorXd ComputeCellForcingTerm(const std::array<Eigen::VectorXd, 3> &forcing_term_values,
                                                  const std::vector<Eigen::MatrixXd> &test_basis_functions_values,
                                                  const Eigen::VectorXd &quadrature_weights) const
    {
        Eigen::MatrixXd rightHandSide =
            test_basis_functions_values[0].transpose() * quadrature_weights.asDiagonal() * forcing_term_values[0];

        for (unsigned int d = 1; d < test_basis_functions_values.size(); d++)
        {
            rightHandSide.noalias() +=
                test_basis_functions_values[d].transpose() * quadrature_weights.asDiagonal() * forcing_term_values[d];
        }

        return rightHandSide;
    }
};
} // namespace Equations
} // namespace PDETools
} // namespace Polydim

#endif
