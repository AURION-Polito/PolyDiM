#include "VEM_Monomials_3D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
//****************************************************************************
VEM_Monomials_Data VEM_Monomials_3D::Compute(const unsigned int polynomial_degree) const
{
    VEM_Monomials_Data data;

    data.Dimension = 3;
    data.PolynomialDegree = polynomial_degree;

    data.NumMonomials = (polynomial_degree + 1) * (polynomial_degree + 2) * (polynomial_degree + 3) / 6;
    data.Exponents.resize(data.NumMonomials);
    data.Exponents[0].setZero(data.Dimension);
    data.DerivativeMatrices.resize(data.Dimension);
    for (unsigned int i = 0; i < data.Dimension; i++)
        data.DerivativeMatrices[i].setZero(data.NumMonomials, data.NumMonomials);
    data.Laplacian.setZero(data.NumMonomials, data.NumMonomials);
    for (unsigned int i = 1; i < data.NumMonomials; i++)
    {
        data.Exponents[i].resize(data.Dimension);
        if (data.Exponents[i - 1](0) == 0)
        {
            if (data.Exponents[i - 1](1) == 0)
                data.Exponents[i] << data.Exponents[i - 1](2) + 1, 0, 0;
            else
                data.Exponents[i] << data.Exponents[i - 1](0) + data.Exponents[i - 1](1) - 1, 0, data.Exponents[i - 1](2) + 1;
        }
        else
            data.Exponents[i] << data.Exponents[i - 1](0) - 1, data.Exponents[i - 1](1) + 1, data.Exponents[i - 1](2);
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
int VEM_Monomials_3D::Index(const VectorXi &exponents) const
{
    int out = -1;
    if (exponents[0] >= 0 && exponents[1] >= 0 && exponents[2] >= 0)
    {
        unsigned int sumExponents = exponents.sum();
        unsigned int sumExponentsIJ = exponents(0) + exponents(1);
        out = (sumExponents + 1) * (sumExponents + 2) * (sumExponents + 3) / 6 -
              (sumExponentsIJ + 1) * (sumExponentsIJ + 2) / 2 + exponents(1);
    }
    return out;
}
//****************************************************************************
std::vector<int> VEM_Monomials_3D::DerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const
{
    std::vector<int> derivativeIndices;
    const VectorXi &expo = data.Exponents[index];
    const unsigned int xyExpoSum = expo(0) + expo(1), xyzExpoSum = xyExpoSum + expo(2);
    const int zDerIndex = index - (xyzExpoSum + 1) * (xyzExpoSum + 2) / 2;
    const int yDerIndex = zDerIndex + xyExpoSum;
    derivativeIndices.resize(3);
    derivativeIndices[0] = (expo(0) > 0) ? yDerIndex + 1 : -1;
    derivativeIndices[1] = (expo(1) > 0) ? yDerIndex : -1;
    derivativeIndices[2] = (expo(2) > 0) ? zDerIndex : -1;
    return derivativeIndices;
}
//****************************************************************************
vector<int> VEM_Monomials_3D::SecondDerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const
{
    vector<int> secondDerivativeIndices;
    const VectorXi &expo = data.Exponents[index];
    const unsigned int xyExpoSum = expo(0) + expo(1), xyzExpoSum = xyExpoSum + expo(2);
    const int zDerIndex = index - (xyzExpoSum + 1) * (xyzExpoSum + 1);
    const int xDerIndex = zDerIndex + 2 * xyExpoSum + 1;
    const int yDerIndex = xDerIndex - 2;
    secondDerivativeIndices.resize(3);
    secondDerivativeIndices[0] = (expo(0) > 1) ? xDerIndex : -1;
    secondDerivativeIndices[1] = (expo(1) > 1) ? yDerIndex : -1;
    secondDerivativeIndices[2] = (expo(2) > 1) ? zDerIndex : -1;
    return secondDerivativeIndices;
}
//****************************************************************************
} // namespace Monomials
} // namespace VEM
} // namespace Polydim
