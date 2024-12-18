#include "VEM_GBasis_3D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
//****************************************************************************
VEM_GBasis_Data VEM_GBasis_3D::Compute(const unsigned int polynomial_degree)
{
    VEM_GBasis_Data data;
    data.Dimension = 3;
    data.PolynomialDegree = polynomial_degree;

    data.monomials_data = monomials.Compute(polynomial_degree + 1);

    data.Nk = (polynomial_degree + 1) * (polynomial_degree + 2) * (polynomial_degree + 3) / 6;
    data.Nkp1 = (polynomial_degree + 2) * (polynomial_degree + 3) * (polynomial_degree + 4) / 6;
    data.Nkm1 = polynomial_degree * (polynomial_degree + 1) * (polynomial_degree + 2) / 6;

    data.DimFirstBasis = data.Nkm1 - (polynomial_degree - 1) * polynomial_degree * (polynomial_degree + 1) / 6;

    data.MapExponents.resize(data.Nkp1);
    for (unsigned int n = 0; n < data.Nkp1; n++)
    {
        data.MapExponents[n].resize(data.Nkp1);
        for (unsigned int m = 0; m < data.Nkp1; m++)
            data.MapExponents[n][m].resize(data.Nkp1);
    }

    for (unsigned int i = 1; i < data.Nkp1; i++)
    {
        if (data.monomials_data.Exponents[i - 1](0) == 0)
        {
            if (data.monomials_data.Exponents[i - 1](1) == 0)
                data.MapExponents[data.monomials_data.Exponents[i - 1](2) + 1][0][0] = i;
            else
                data.MapExponents[data.monomials_data.Exponents[i - 1](0) + data.monomials_data.Exponents[i - 1](1) - 1]
                                 [0][data.monomials_data.Exponents[i - 1](2) + 1] = i;
        }
        else
            data.MapExponents[data.monomials_data.Exponents[i - 1](0) - 1][data.monomials_data.Exponents[i - 1](1) + 1]
                             [data.monomials_data.Exponents[i - 1](2)] = i;
    }

    data.VectorDecomposition.resize(data.Dimension);

    for (unsigned int i = 0; i < data.Dimension; i++)
    {
        data.VectorDecomposition[i].resize(4);
        data.VectorDecomposition[i][0].setZero(data.Nk, data.Nkp1);
        data.VectorDecomposition[i][1].setZero(data.Nk, data.DimFirstBasis);
        data.VectorDecomposition[i][2].setZero(data.Nk, data.Nkm1);
        data.VectorDecomposition[i][3].setZero(data.Nk, data.Nkm1);
    }

    data.MapFirstGroupVectorDecomposition.resize(data.Nkm1, -1);
    for (unsigned int j = 0; j < data.Nkm1; j++)
    {
        const VectorXi &expo = data.monomials_data.Exponents[j];
        if (expo(0) == 0)
        {
            data.MapFirstGroupVectorDecomposition[j] = data.SaveFirstGroupVectorDecomposition.size();
            data.SaveFirstGroupVectorDecomposition.push_back(j);
        }
    }

    for (unsigned int j = 0; j < data.Nk; j++)
    {
        const VectorXi &expo = data.monomials_data.Exponents[j];
        const vector<Vector4i> VectorDecompositionIndex = VectorDecompositionIndices(data, expo);
        const double invSumExponent = 1.0 / (data.monomials_data.Exponents[j].sum() + 1.0);

        data.VectorDecomposition[0][0](j, VectorDecompositionIndex[0](0)) = invSumExponent;
        if (expo(2) > 0)
            data.VectorDecomposition[0][2](j, VectorDecompositionIndex[0](2)) = -((double)expo(2)) * invSumExponent;
        if (expo(1) > 0)
            data.VectorDecomposition[0][3](j, VectorDecompositionIndex[0](3)) = ((double)expo(1)) * invSumExponent;

        data.VectorDecomposition[1][0](j, VectorDecompositionIndex[1](0)) = invSumExponent;
        if (expo(2) > 0 && expo(0) == 0)
            data.VectorDecomposition[1][1](j, data.MapFirstGroupVectorDecomposition[VectorDecompositionIndex[1](1)]) =
                ((double)expo(2)) * invSumExponent;
        if (expo(0) > 0 && expo(2) > 0)
            data.VectorDecomposition[1][2](j, VectorDecompositionIndex[1](2)) = -((double)expo(2)) * invSumExponent;
        if (expo(0) > 0)
            data.VectorDecomposition[1][3](j, VectorDecompositionIndex[1](3)) =
                -((double)(expo(0) + expo(2))) * invSumExponent;

        data.VectorDecomposition[2][0](j, VectorDecompositionIndex[2](0)) = invSumExponent;
        if (expo(1) > 0 && expo(0) == 0)
            data.VectorDecomposition[2][1](j, data.MapFirstGroupVectorDecomposition[VectorDecompositionIndex[2](1)]) =
                -((double)expo(1)) * invSumExponent;
        if (expo(0) > 0)
            data.VectorDecomposition[2][2](j, VectorDecompositionIndex[2](2)) =
                ((double)(expo(0) + expo(1))) * invSumExponent;
        if (expo(0) > 0 && expo(1) > 0)
            data.VectorDecomposition[2][3](j, VectorDecompositionIndex[2](3)) = ((double)expo(1)) * invSumExponent;
    }

    return data;
}
//****************************************************************************
vector<Vector4i> VEM_GBasis_3D::VectorDecompositionIndices(const VEM_GBasis_Data &data, const VectorXi &expo) const
{
    vector<Vector4i> vectorDecompositionIndices(3, Vector4i::Constant(-1));

    vectorDecompositionIndices[0](0) = data.MapExponents[expo(0) + 1][expo(1)][expo(2)];
    vectorDecompositionIndices[0](2) = (expo(2) > 0) ? data.MapExponents[expo(0)][expo(1)][expo(2) - 1] : -1;
    vectorDecompositionIndices[0](3) = (expo(1) > 0) ? data.MapExponents[expo(0)][expo(1) - 1][expo(2)] : -1;

    vectorDecompositionIndices[1](0) = data.MapExponents[expo(0)][expo(1) + 1][expo(2)];
    vectorDecompositionIndices[1](1) =
        (expo(0) == 0 && expo(2) > 0) ? data.MapExponents[expo(0)][expo(1)][expo(2) - 1] : -1;
    vectorDecompositionIndices[1](2) =
        (expo(0) > 0 && expo(2) > 0) ? data.MapExponents[expo(0) - 1][expo(1) + 1][expo(2) - 1] : -1;
    vectorDecompositionIndices[1](3) = (expo(0) > 0) ? data.MapExponents[expo(0) - 1][expo(1)][expo(2)] : -1;

    vectorDecompositionIndices[2](0) = data.MapExponents[expo(0)][expo(1)][expo(2) + 1];
    vectorDecompositionIndices[2](1) =
        (expo(0) == 0 && expo(1) > 0) ? data.MapExponents[expo(0)][expo(1) - 1][expo(2)] : -1;
    vectorDecompositionIndices[2](2) = (expo(0) > 0) ? data.MapExponents[expo(0) - 1][expo(1)][expo(2)] : -1;
    vectorDecompositionIndices[2](3) =
        (expo(0) > 0 && expo(1) > 0) ? data.MapExponents[expo(0) - 1][expo(1) - 1][expo(2) + 1] : -1;

    return vectorDecompositionIndices;
}
//****************************************************************************
vector<MatrixXd> VEM_GBasis_3D::VanderGBigOPlus(const VEM_GBasis_Data &data, const MatrixXd &vander) const
{

    const unsigned int numPoints = vander.rows();
    vector<MatrixXd> vanderBigOPlus(3, MatrixXd::Zero(numPoints, data.DimFirstBasis + 2 * data.Nkm1));

    if (vander.cols() < 4)
        return vanderBigOPlus;

    const MatrixXd tmpZ = vander.leftCols(data.Nkm1).array().colwise() * vander.col(3).array();
    const MatrixXd tmpY = vander.leftCols(data.Nkm1).array().colwise() * vander.col(2).array();
    const MatrixXd tmpX = vander.leftCols(data.Nkm1).array().colwise() * vander.col(1).array();

    unsigned int count = 0;
    for (unsigned int i : data.SaveFirstGroupVectorDecomposition)
    {
        vanderBigOPlus[1].col(count) = tmpZ.col(i);
        vanderBigOPlus[2].col(count++) = -tmpY.col(i);
    }

    vanderBigOPlus[0].middleCols(data.DimFirstBasis, data.Nkm1) = -tmpZ;
    vanderBigOPlus[1].middleCols(data.DimFirstBasis, data.Nkm1) = MatrixXd::Zero(numPoints, data.Nkm1);
    vanderBigOPlus[2].middleCols(data.DimFirstBasis, data.Nkm1) = tmpX;

    vanderBigOPlus[0].middleCols(data.DimFirstBasis + data.Nkm1, data.Nkm1) = tmpY;
    vanderBigOPlus[1].middleCols(data.DimFirstBasis + data.Nkm1, data.Nkm1) = -tmpX;
    vanderBigOPlus[2].middleCols(data.DimFirstBasis + data.Nkm1, data.Nkm1) = MatrixXd::Zero(numPoints, data.Nkm1);

    return vanderBigOPlus;
}
//****************************************************************************
} // namespace Monomials
} // namespace VEM
} // namespace Polydim
