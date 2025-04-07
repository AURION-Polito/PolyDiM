#include "VEM_GBasis_2D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
//****************************************************************************
VEM_GBasis_Data VEM_GBasis_2D::Compute(const unsigned int polynomial_degree)
{
    VEM_GBasis_Data data;
    data.Dimension = 2;
    data.PolynomialDegree = polynomial_degree;

    data.monomials_data = monomials.Compute(polynomial_degree + 1);

    data.Nk = (polynomial_degree + 1) * (polynomial_degree + 2) * 0.5;
    data.Nkp1 = (polynomial_degree + 2) * (polynomial_degree + 3) * 0.5;
    data.Nkm1 = polynomial_degree * (polynomial_degree + 1) * 0.5;

    data.MatrixExponents = MatrixXi::Constant(data.Nkp1, data.Nkp1, -1);
    data.MatrixExponents(0, 0) = 0;
    for (unsigned int i = 1; i < data.Nkp1; i++)
    {
        if (data.monomials_data.Exponents[i - 1](0) == 0)
            data.MatrixExponents(data.monomials_data.Exponents[i - 1](1) + 1, 0) = i;
        else
            data.MatrixExponents(data.monomials_data.Exponents[i - 1](0) - 1, data.monomials_data.Exponents[i - 1](1) + 1) = i;
    }

    data.VectorDecomposition.resize(data.Dimension);

    for (unsigned int i = 0; i < data.Dimension; i++)
    {
        data.VectorDecomposition[i].resize(data.Dimension);
        data.VectorDecomposition[i][0].setZero(data.Nk, data.Nkp1);
        data.VectorDecomposition[i][1].setZero(data.Nk, data.Nkm1);
    }

    for (unsigned int j = 0; j < data.Nk; j++)
    {
        const VectorXi &expo = data.monomials_data.Exponents[j];
        const vector<Vector2i> VectorDecompositionIndex = VectorDecompositionIndices(data, expo);

        data.VectorDecomposition[0][0](j, VectorDecompositionIndex[0](0)) = 1.0 / (data.monomials_data.Exponents[j].sum() + 1);
        if (expo(1) > 0)
            data.VectorDecomposition[0][1](j, VectorDecompositionIndex[0](1)) =
                ((double)expo(1)) / (data.monomials_data.Exponents[j].sum() + 1);

        data.VectorDecomposition[1][0](j, VectorDecompositionIndex[1](0)) = 1.0 / (data.monomials_data.Exponents[j].sum() + 1);
        if (expo(0) > 0)
            data.VectorDecomposition[1][1](j, VectorDecompositionIndex[1](1)) =
                -((double)expo(0)) / (data.monomials_data.Exponents[j].sum() + 1);
    }

    return data;
}
//****************************************************************************
vector<Vector2i> VEM_GBasis_2D::VectorDecompositionIndices(const VEM_GBasis_Data &data, const VectorXi &expo) const
{
    vector<Vector2i> vectorDecompositionIndices(2);

    vectorDecompositionIndices[0](0) = data.MatrixExponents(expo(0) + 1, expo(1));
    vectorDecompositionIndices[0](1) = (expo(1) > 0) ? data.MatrixExponents(expo(0), expo(1) - 1) : -1;

    vectorDecompositionIndices[1](0) = data.MatrixExponents(expo(0), expo(1) + 1);
    vectorDecompositionIndices[1](1) = (expo(0) > 0) ? data.MatrixExponents(expo(0) - 1, expo(1)) : -1;

    return vectorDecompositionIndices;
}
//****************************************************************************
} // namespace Utilities
} // namespace VEM
} // namespace Polydim
