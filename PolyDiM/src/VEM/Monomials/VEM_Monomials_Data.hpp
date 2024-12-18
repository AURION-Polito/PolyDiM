#ifndef __VEM_Monomials_Data_HPP
#define __VEM_Monomials_Data_HPP

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
struct VEM_Monomials_Data
{
    unsigned int PolynomialDegree;                   ///< Monomial space order
    unsigned int Dimension;                          ///< The geometric dimension
    unsigned int NumMonomials;                       ///< Number of monomials in the basis.
    std::vector<Eigen::VectorXi> Exponents;          ///< Table of exponents of each monomial.
    std::vector<Eigen::MatrixXd> DerivativeMatrices; ///< Matrices used to compute derivatives of monomials.
    Eigen::MatrixXd Laplacian;                       ///< Matrix used to compute the laplacian of monomials.
};
} // namespace Monomials
} // namespace VEM
} // namespace Polydim

#endif
