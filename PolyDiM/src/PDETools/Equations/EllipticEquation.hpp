#ifndef __EllipticEquation_H
#define __EllipticEquation_H

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
          Eigen::MatrixXd ComputeCellDiffusionMatrix(const Eigen::VectorXd& diffusion_term_values,
                                                     const std::vector<Eigen::MatrixXd>& basis_functions_derivative_values,
                                                     const Eigen::VectorXd& quadrature_weights) const
          {
            Eigen::MatrixXd cell_matrix =
                basis_functions_derivative_values.at(0).transpose() *
                quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                basis_functions_derivative_values.at(0);

            for (unsigned int d = 1; d < basis_functions_derivative_values.size(); ++d)
            {
              cell_matrix.noalias() +=
                  basis_functions_derivative_values.at(1).transpose() *
                  quadrature_weights.cwiseProduct(diffusion_term_values).asDiagonal() *
                  basis_functions_derivative_values.at(1);
            }

            return cell_matrix;
          }

          inline Eigen::MatrixXd ComputeCellReactionMatrix(const Eigen::VectorXd& reaction_term_values,
                                                           const Eigen::MatrixXd& basis_functions_values,
                                                           const Eigen::VectorXd& quadrature_weights) const
          {

            return basis_functions_values.transpose() *
                quadrature_weights.cwiseProduct(reaction_term_values).asDiagonal() *
                basis_functions_values;
          }

          Eigen::MatrixXd ComputeCellAdvectionMatrix(const std::vector<Eigen::VectorXd>& advection_term_values,
                                                     const Eigen::MatrixXd& basis_functions_values,
                                                     const std::vector<Eigen::MatrixXd>& basis_functions_derivative_values,
                                                     const Eigen::VectorXd& quadrature_weights) const
          {
            Eigen::MatrixXd cell_matrix = basis_functions_values.transpose() *
                                          quadrature_weights.cwiseProduct(advection_term_values[0]).asDiagonal() *
                basis_functions_derivative_values[0];

            for(unsigned int d = 1; d < basis_functions_derivative_values.size(); ++d)
            {
              cell_matrix.noalias() +=
                  basis_functions_values.transpose() *
                  quadrature_weights.cwiseProduct(advection_term_values.at(d)).asDiagonal() *
                  basis_functions_derivative_values.at(d);
            }

            return cell_matrix;
          }

          inline Eigen::VectorXd ComputeCellForcingTerm(const Eigen::VectorXd& forcingTermValues,
                                                        const Eigen::MatrixXd& basis_functions_values,
                                                        const Eigen::VectorXd& quadrature_weights) const
          {
            return basis_functions_values.transpose() *
                quadrature_weights.asDiagonal() *
                forcingTermValues;
          }
      };
    }
  }
}

#endif
