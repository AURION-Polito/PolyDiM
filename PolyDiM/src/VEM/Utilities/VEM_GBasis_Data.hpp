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

#ifndef __VEM_GBasis_Data_HPP
#define __VEM_GBasis_Data_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
struct VEM_GBasis_Data
{
    unsigned int PolynomialDegree;
    unsigned int Dimension;

    VEM_Monomials_Data monomials_data;

    unsigned int Nk;
    unsigned int Nkm1;
    unsigned int Nkp1;
    unsigned int NkGBigOPlus;
    unsigned int NkGNabla;
    std::vector<std::vector<Eigen::MatrixXd>> VectorDecomposition;

    Eigen::MatrixXi MatrixExponents;

    unsigned int DimFirstBasis;
    std::vector<std::vector<std::vector<unsigned int>>> MapExponents;
    std::vector<unsigned int> SaveFirstGroupVectorDecomposition;
    std::vector<int> MapFirstGroupVectorDecomposition;
};
} // namespace Utilities
} // namespace VEM
} // namespace Polydim

#endif
