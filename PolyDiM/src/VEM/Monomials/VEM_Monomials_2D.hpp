#ifndef __VEM_Monomials_VEM_Monomials_2D_HPP
#define __VEM_Monomials_VEM_Monomials_2D_HPP

#include "VEM_Monomials_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
class VEM_Monomials_2D final
{
  private:
    VEM_Monomials_Utilities<2> utilities;

  public:
    VEM_Monomials_Data Compute(const unsigned int polynomial_degree) const;

    inline Eigen::MatrixXi Exponents(const VEM_Monomials_Data &data) const
    {
        return utilities.Exponents(data);
    }

    /// \param The required derivative.
    /// \returns A const reference to \ref derivativeMatrices[i].
    inline Eigen::MatrixXd DerivativeMatrix(const VEM_Monomials_Data &data, const unsigned int &i) const
    {
        return data.DerivativeMatrices[i];
    }

    /// \returns A const reference to the x-derivative matrix.
    inline Eigen::MatrixXd D_x(const VEM_Monomials_Data &data) const
    {
        return data.DerivativeMatrices[0];
    }

    /// \returns A const reference to the y-derivative matrix.
    inline Eigen::MatrixXd D_y(const VEM_Monomials_Data &data) const
    {
        return data.DerivativeMatrices[1];
    }

    int Index(const Eigen::VectorXi &exponents) const;

    std::vector<int> DerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const;

    std::vector<int> SecondDerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const;

    inline Eigen::MatrixXd Vander(const VEM_Monomials_Data &data,
                                  const std::vector<Eigen::VectorXd> &points,
                                  const Eigen::VectorXd &centroid,
                                  const double &diam) const
    {
        return utilities.Vander(data, points, centroid, diam);
    }
    inline Eigen::MatrixXd Vander(const VEM_Monomials_Data &data,
                                  const Eigen::MatrixXd &points,
                                  const Eigen::Vector3d &centroid,
                                  const double &diam) const
    {
        return utilities.Vander(data, points, centroid, diam);
    }
    inline std::vector<Eigen::MatrixXd> VanderDerivatives(const VEM_Monomials_Data &data,
                                                          const Eigen::MatrixXd &vander,
                                                          const double &diam) const
    {
        return utilities.VanderDerivatives(data, (*this), vander, diam);
    }
    inline Eigen::MatrixXd VanderLaplacian(const VEM_Monomials_Data &data,
                                           const Eigen::MatrixXd &vander,
                                           const double &diam) const
    {
        return utilities.VanderLaplacian(data, (*this), vander, diam);
    }
};
} // namespace Monomials
} // namespace VEM
} // namespace Polydim

#endif
