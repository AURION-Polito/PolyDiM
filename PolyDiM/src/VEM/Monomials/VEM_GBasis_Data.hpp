#ifndef __VEM_GBasis_Data_HPP
#define __VEM_GBasis_Data_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
struct VEM_GBasis_Data
{
    unsigned int PolynomialDegree; ///< Monomial space order
    unsigned int Dimension;        ///< The geometric dimension

    VEM_Monomials_Data monomials_data;

    unsigned int Nk;   ///< Number of monomials in the basis.
    unsigned int Nkm1; ///< Number of monomials in the basis.
    unsigned int Nkp1; ///< Number of monomials in the basis.
    unsigned int NkGBigOPlus;
    unsigned int NkGNabla;
    std::vector<std::vector<Eigen::MatrixXd>> VectorDecomposition; ///< Matrix used to compute the laplacian of
                                                                   ///< monomials.

    Eigen::MatrixXi MatrixExponents; ///< Table of exponents of each monomial.

    unsigned int DimFirstBasis;
    std::vector<std::vector<std::vector<unsigned int>>> MapExponents;
    std::vector<unsigned int> SaveFirstGroupVectorDecomposition;
    std::vector<int> MapFirstGroupVectorDecomposition;
};
} // namespace Monomials
} // namespace VEM
} // namespace Polydim

#endif
