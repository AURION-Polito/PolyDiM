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
    /// @brief Local diffusion matrix with a scalar diffusion coefficient (Petrov–Galerkin).
    ///
    /// Computes \f$\int_E \mu\, \nabla u \cdot \nabla v\f$ by summing the products of the
    /// trial and test derivatives over the spatial directions, weighted by the scalar
    /// diffusion field \f$\mu\f$ and the quadrature weights.
    ///
    /// @param diffusion_term_values                    Scalar diffusion coefficient \f$\mu\f$ at the quadrature points.
    /// @param trial_basis_functions_derivative_values  Trial-space derivatives, one matrix per spatial direction.
    /// @param test_basis_functions_derivative_values   Test-space derivatives, one matrix per spatial direction.
    /// @param quadrature_weights                       Quadrature weights.
    /// @return The local diffusion matrix (test DOFs \f$\times\f$ trial DOFs).
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

    /// @brief Local diffusion matrix with a scalar diffusion coefficient (Galerkin).
    ///
    /// Galerkin specialization of the Petrov–Galerkin overload, using the same basis
    /// for trial and test spaces: \f$\int_E \mu\, \nabla u \cdot \nabla v\f$.
    ///
    /// @param diffusion_term_values           Scalar diffusion coefficient \f$\mu\f$ at the quadrature points.
    /// @param basis_functions_derivative_values Basis-function derivatives, one matrix per spatial direction.
    /// @param quadrature_weights              Quadrature weights.
    /// @return The local diffusion matrix.
    Eigen::MatrixXd ComputeCellDiffusionMatrix(const Eigen::VectorXd &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        return ComputeCellDiffusionMatrix(diffusion_term_values, basis_functions_derivative_values, basis_functions_derivative_values, quadrature_weights);
    }

    /// @brief Local diffusion matrix with a full diffusion tensor (Petrov–Galerkin).
    ///
    /// Computes \f$\int_E (\mathbf{K}\, \nabla u) \cdot \nabla v\f$ for a full
    /// (up to \f$3\times 3\f$) diffusion tensor \f$\mathbf{K}\f$, summing over both
    /// derivative directions. The tensor is passed in column-major order, i.e. the
    /// entry \f$K_{d_1 d_2}\f$ is stored at index \f$d_1 + 3\,d_2\f$.
    ///
    /// @param diffusion_term_values                    Diffusion tensor entries \f$K_{d_1 d_2}\f$ (column-major, length
    /// 9) at the quadrature points.
    /// @param trial_basis_functions_derivative_values  Trial-space derivatives, one matrix per spatial direction.
    /// @param test_basis_functions_derivative_values   Test-space derivatives, one matrix per spatial direction.
    /// @param quadrature_weights                       Quadrature weights.
    /// @return The local diffusion matrix (test DOFs \f$\times\f$ trial DOFs).
    Eigen::MatrixXd ComputeCellDiffusionMatrix(const std::array<Eigen::VectorXd, 9> &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &trial_basis_functions_derivative_values,
                                               const std::vector<Eigen::MatrixXd> &test_basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        const unsigned int dimension = trial_basis_functions_derivative_values.size();

        Eigen::MatrixXd cell_matrix = Eigen::MatrixXd::Zero(test_basis_functions_derivative_values.at(0).cols(),
                                                            trial_basis_functions_derivative_values.at(0).cols());
        for (unsigned int d1 = 0; d1 < dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < dimension; d2++)
            {
                cell_matrix.noalias() += test_basis_functions_derivative_values.at(d1).transpose() *
                                         quadrature_weights.cwiseProduct(diffusion_term_values.at(d1 + 3 * d2)).asDiagonal() *
                                         trial_basis_functions_derivative_values.at(d2);
            }
        }
        return cell_matrix;
    }

    /// @brief Local diffusion matrix with a full diffusion tensor (Galerkin).
    ///
    /// Galerkin specialization of the tensor Petrov–Galerkin overload:
    /// \f$\int_E (\mathbf{K}\, \nabla u) \cdot \nabla v\f$ with a shared basis.
    ///
    /// @param diffusion_term_values           Diffusion tensor entries (column-major, length 9) at the quadrature
    /// points.
    /// @param basis_functions_derivative_values Basis-function derivatives, one matrix per spatial direction.
    /// @param quadrature_weights              Quadrature weights.
    /// @return The local diffusion matrix.
    Eigen::MatrixXd ComputeCellDiffusionMatrix(const std::array<Eigen::VectorXd, 9> &diffusion_term_values,
                                               const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                               const Eigen::VectorXd &quadrature_weights) const
    {
        return ComputeCellDiffusionMatrix(diffusion_term_values, basis_functions_derivative_values, basis_functions_derivative_values, quadrature_weights);
    }

    /// @brief Local reaction (mass-like) matrix (Galerkin).
    ///
    /// Computes \f$\int_E \sigma\, u\, v\f$ with reaction coefficient \f$\sigma\f$,
    /// using the same basis for trial and test spaces.
    ///
    /// @param reaction_term_values   Reaction coefficient \f$\sigma\f$ at the quadrature points.
    /// @param basis_functions_values Basis-function values at the quadrature points.
    /// @param quadrature_weights     Quadrature weights.
    /// @return The local reaction matrix.
    inline Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd &reaction_term_values,
                                                     const Eigen::MatrixXd &basis_functions_values,
                                                     const Eigen::VectorXd &quadrature_weights) const
    {
        return basis_functions_values.transpose() * quadrature_weights.cwiseProduct(reaction_term_values).asDiagonal() * basis_functions_values;
    }

    /// @brief Local reaction (mass-like) matrix (Petrov–Galerkin).
    ///
    /// Computes \f$\int_E \sigma\, u\, v\f$ with distinct trial and test spaces.
    ///
    /// @param reaction_term_values         Reaction coefficient \f$\sigma\f$ at the quadrature points.
    /// @param trial_basis_functions_values Trial-space values at the quadrature points.
    /// @param test_basis_functions_values  Test-space values at the quadrature points.
    /// @param quadrature_weights           Quadrature weights.
    /// @return The local reaction matrix (test DOFs \f$\times\f$ trial DOFs).
    inline Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd &reaction_term_values,
                                                     const Eigen::MatrixXd &trial_basis_functions_values,
                                                     const Eigen::MatrixXd &test_basis_functions_values,
                                                     const Eigen::VectorXd &quadrature_weights) const
    {
        return test_basis_functions_values.transpose() *
               quadrature_weights.cwiseProduct(reaction_term_values).asDiagonal() * trial_basis_functions_values;
    }

    /// @brief Local advection matrix (Petrov–Galerkin).
    ///
    /// Computes \f$\int_E (\boldsymbol{\beta} \cdot \nabla u)\, v\f$ for an advection
    /// field \f$\boldsymbol{\beta}\f$ (up to 3 components), pairing the test-function
    /// values with the trial-function derivatives summed over the spatial directions.
    ///
    /// @param advection_term_values                    Advection field components \f$\beta_d\f$ at the quadrature
    /// points.
    /// @param test_basis_functions_values              Test-space values at the quadrature points.
    /// @param trial_basis_functions_derivative_values  Trial-space derivatives, one matrix per spatial direction.
    /// @param quadrature_weights                       Quadrature weights.
    /// @return The local advection matrix (test DOFs \f$\times\f$ trial DOFs).
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

    /// @brief Local forcing term for a scalar problem.
    ///
    /// Computes \f$\int_E f\, v\f$ with scalar source \f$f\f$.
    ///
    /// @param forcing_term_values         Source term \f$f\f$ at the quadrature points.
    /// @param test_basis_functions_values Test-space values at the quadrature points.
    /// @param quadrature_weights          Quadrature weights.
    /// @return The local right-hand-side vector.
    inline Eigen::VectorXd ComputeCellForcingTerm(const Eigen::VectorXd &forcing_term_values,
                                                  const Eigen::MatrixXd &test_basis_functions_values,
                                                  const Eigen::VectorXd &quadrature_weights) const
    {
        return test_basis_functions_values.transpose() * quadrature_weights.asDiagonal() * forcing_term_values;
    }

    /// @brief Local forcing term for a vector-valued problem.
    ///
    /// Computes \f$\int_E \boldsymbol{f} \cdot \boldsymbol{v}\f$ for a vector source
    /// \f$\boldsymbol{f}\f$ (up to 3 components), summing the contribution of each
    /// component with its corresponding test-function component.
    ///
    /// @param forcing_term_values         Source components \f$f_d\f$ at the quadrature points.
    /// @param test_basis_functions_values Test-space values per component, one matrix per direction.
    /// @param quadrature_weights          Quadrature weights.
    /// @return The local right-hand-side vector.
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
