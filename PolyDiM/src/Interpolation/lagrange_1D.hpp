#ifndef __Interpolation_Lagrange_1D_H
#define __Interpolation_Lagrange_1D_H

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
  namespace Interpolation
  {
    namespace Lagrange
    {
      Eigen::VectorXd Lagrange_1D_cofficients(const Eigen::VectorXd& interpolation_points_x);
      Eigen::MatrixXd Lagrange_1D_values(const Eigen::VectorXd& interpolation_points_x,
                                         const Eigen::VectorXd& lagrange_1D_coefficients,
                                         const Eigen::VectorXd& evaluation_points_x);
      Eigen::MatrixXd Lagrange_1D_derivative_values(const Eigen::VectorXd& interpolation_points_x,
                                                    const Eigen::VectorXd& lagrange_1D_coefficients,
                                                    const Eigen::VectorXd& evaluation_points_x);


    } // namespace Lagrange
  } // namespace Interpolation
} // namespace Polydim

#endif
