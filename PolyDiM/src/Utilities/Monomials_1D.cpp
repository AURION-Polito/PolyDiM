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

#include "Monomials_1D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace Utilities
{
//****************************************************************************
Monomials_Data Monomials_1D::Compute(const unsigned int polynomial_degree) const
{
    Monomials_Data data;

    data.Dimension = 1;
    data.PolynomialDegree = polynomial_degree;
    data.NumMonomials = polynomial_degree + 1;
    data.Exponents.resize(data.NumMonomials);
    data.Exponents[0].setZero(data.Dimension);
    data.DerivativeMatrices.resize(data.Dimension);
    for (unsigned int i = 0; i < data.Dimension; i++)
        data.DerivativeMatrices[i].setZero(data.NumMonomials, data.NumMonomials);
    data.Laplacian.setZero(data.NumMonomials, data.NumMonomials);
    for (unsigned int i = 1; i < data.NumMonomials; i++)
    {
        data.Exponents[i].resize(data.Dimension);
        data.Exponents[i] << data.Exponents[i - 1][0] + 1;
        vector<int> derIndices = DerivativeIndices(data, i);
        vector<int> secondDerIndices = SecondDerivativeIndices(data, i);
        for (unsigned int j = 0; j < data.Dimension; ++j)
        {
            if (derIndices[j] >= 0)
                data.DerivativeMatrices[j](i, derIndices[j]) = data.Exponents[i](j);
            if (secondDerIndices[j] >= 0)
                data.Laplacian(i, secondDerIndices[j]) = data.Exponents[i](j) * (data.Exponents[i](j) - 1);
        }
    }

    return data;
}
//****************************************************************************
int Monomials_1D::Index(const VectorXi &exponents) const
{
    int out = -1;
    if (exponents[0] >= 0)
        out = exponents.sum();
    return out;
}
//****************************************************************************
vector<int> Monomials_1D::DerivativeIndices(const Monomials_Data &data, const unsigned int &index) const
{
    vector<int> derivativeIndices;

    const VectorXi &expo = data.Exponents[index];
    const int xDerIndex = index - 1;
    derivativeIndices.resize(1);
    derivativeIndices[0] = (expo(0) > 0) ? xDerIndex : -1;

    return derivativeIndices;
}
//****************************************************************************
vector<int> Monomials_1D::SecondDerivativeIndices(const Monomials_Data &data, const unsigned int &index) const
{
    vector<int> secondDerivativeIndices;
    const VectorXi &expo = data.Exponents[index];
    const int xxDerIndex = index - 2;

    secondDerivativeIndices.resize(1);
    secondDerivativeIndices[0] = (expo(0) > 1) ? xxDerIndex : -1;
    return secondDerivativeIndices;
}
//****************************************************************************
} // namespace Utilities
} // namespace Polydim
