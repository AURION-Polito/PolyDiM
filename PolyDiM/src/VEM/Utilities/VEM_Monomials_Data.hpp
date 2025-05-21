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

#ifndef __VEM_Monomials_Data_HPP
#define __VEM_Monomials_Data_HPP

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace Utilities
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
} // namespace Utilities
} // namespace VEM
} // namespace Polydim

#endif
