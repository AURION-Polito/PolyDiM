#ifndef __VEM_Monomials_VEM_Monomials_1D_HPP
#define __VEM_Monomials_VEM_Monomials_1D_HPP

#include "VEM_Monomials_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
class VEM_Monomials_1D final
{
  private:
    VEM_Monomials_Utilities<1> utilities;

  public:
    /// \brief Define the multi-indices \f$\alpha\f$ data for \f$(x,y,z)\f$ and the other coefficient matrices
    /// to build efficiently the monomial Vandermonde matrix and
    /// the Vandermonde Matrix for both monomials derivatives and laplacian.
    /// \param polynomial_degree: The polynomial degree required.
    /// \returns A struct of type \ref VEM_Monomials_Data
    VEM_Monomials_Data Compute(const unsigned int polynomial_degree) const;

    /// \param A const reference to an object of type \ref VEM_Monomials_Data.
    /// \returns The an Eigen::MatrixXi containing the \ref VEM_Monomials_Data::Exponents
    inline Eigen::MatrixXi Exponents(const VEM_Monomials_Data &data) const
    {
        return utilities.Exponents(data);
    }

    /// \param A const reference to an object of type \ref VEM_Monomials_Data.
    /// \param An unsigned int that expresses the required derivative:
    /// \f$0\f$ for the \f$x\f$-derivative
    /// \returns The \f$i\f$-th entry of \ref VEM_Monomials_Data::DerivativeMatrices.
    inline Eigen::MatrixXd DerivativeMatrix(const VEM_Monomials_Data &data, const unsigned int &i) const
    {
        return data.DerivativeMatrices[i];
    }

    /// \param A const reference to an object of type \ref VEM_Monomials_Data.
    /// \returns A const reference to the x-derivative matrix.
    inline Eigen::MatrixXd D_x(const VEM_Monomials_Data &data) const
    {
        return data.DerivativeMatrices[0];
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
