#include "VEM_MCC_Utilities.hpp"

#include "CommonUtilities.hpp"

using namespace Eigen;
using namespace std;

namespace Polydim
{
namespace VEM
{
namespace MCC
{
template struct VEM_MCC_Utilities<2>;
template struct VEM_MCC_Utilities<3>;
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_MCC_Utilities<dimension>::ComputeDofiDofiStabilizationMatrix(const MatrixXd &projector,
                                                                          const double &coefficient,
                                                                          const Eigen::MatrixXd &DMatrix) const
{
    MatrixXd stabMatrixPi0k = DMatrix * projector;
    stabMatrixPi0k.diagonal().array() -= 1;
    // stabMatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
    stabMatrixPi0k = coefficient * stabMatrixPi0k.transpose() * stabMatrixPi0k;

    return stabMatrixPi0k;
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_MCC_Utilities<dimension>::ComputePolynomialBasisDofs(const double &polytopeMeasure,
                                                                  const unsigned int &order,
                                                                  const unsigned int &Nk,
                                                                  const unsigned int &NumBoundaryBasisFunctions,
                                                                  const unsigned int &NumNablaInternalBasisFunctions,
                                                                  const unsigned int &NumBigOPlusInternalBasisFunctions,
                                                                  const unsigned int &NumBasisFunctions,
                                                                  const Eigen::MatrixXd &GkVanderBoundaryTimesNormal,
                                                                  const Eigen::MatrixXd &Gmatrix) const
{
    MatrixXd polynomialBasisDofs = MatrixXd::Zero(NumBasisFunctions, dimension * Nk);

    polynomialBasisDofs.topRows(NumBoundaryBasisFunctions) = GkVanderBoundaryTimesNormal;

    if (order > 0)
    {
        polynomialBasisDofs.block(NumBoundaryBasisFunctions, 0, NumNablaInternalBasisFunctions, dimension * Nk) =
            (1.0 / polytopeMeasure) * Gmatrix.topRows(NumNablaInternalBasisFunctions);

        polynomialBasisDofs.bottomRows(NumBigOPlusInternalBasisFunctions) =
            (1.0 / polytopeMeasure) * Gmatrix.bottomRows(NumBigOPlusInternalBasisFunctions);
    }
    return polynomialBasisDofs;
}
//****************************************************************************
template <unsigned short dimension>
std::vector<MatrixXd> VEM_MCC_Utilities<dimension>::ComputeBasisFunctionsValues(const Eigen::MatrixXd &projector,
                                                                                const Eigen::MatrixXd &GVander) const
{
    const unsigned int numQuadrature = GVander.cols() / dimension;
    const Eigen::MatrixXd temp = GVander.transpose() * projector;
    std::vector<Eigen::MatrixXd> result(dimension, Eigen::MatrixXd::Zero(dimension, projector.cols()));

    for (unsigned int d = 0; d < dimension; d++)
        result[d] = temp.middleRows(numQuadrature * d, numQuadrature);

    return result;
}
//****************************************************************************
template <unsigned short dimension>
void VEM_MCC_Utilities<dimension>::MonomialTraceOnEdges(const unsigned int &polynomialDegree,
                                                        const Eigen::MatrixXd &polygonVertices,
                                                        const double &polygonDiameter,
                                                        const Eigen::Vector3d &polygonCentroid,
                                                        const std::vector<bool> &edgeDirections,
                                                        const Eigen::MatrixXd &edgeTangents,
                                                        std::vector<MatrixXd> &Cmatrixkp1) const
{
    std::vector<Eigen::VectorXd> BinomialCoefficients(polynomialDegree + 2);

    for (unsigned int n = 0; n <= (polynomialDegree + 1); n++)
    {
        BinomialCoefficients[n].resize(n + 2);
        for (unsigned int k = 0; k <= n; k++)
        {
            BinomialCoefficients[n][k] = Gedim::Utilities::BinomialCoefficient(n, k);
        }
    }

    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    const unsigned int order_1D = (polynomialDegree + 2);
    const unsigned int nkp1 = (polynomialDegree + 2) * (polynomialDegree + 3) / 2;

    Cmatrixkp1.resize(numEdges);

    for (unsigned int e = 0; e < numEdges; ++e)
    {
        Cmatrixkp1[e] = MatrixXd::Zero(nkp1, order_1D);

        const Vector3d &edgeStart = edgeDirections[e] ? polygonVertices.col(e) : polygonVertices.col((e + 1) % numVertices);
        const Vector3d &edgeTangent = edgeTangents.col(e);
        const double direction = edgeDirections[e] ? 1.0 : -1.0;

        double invDiameter = 1.0 / polygonDiameter;
        double XAminusXC = edgeStart(0) - polygonCentroid(0);
        double YAminusYC = edgeStart(1) - polygonCentroid(1);
        double X = XAminusXC * invDiameter;
        double Y = YAminusYC * invDiameter;
        double diffX = direction * edgeTangent(0) / XAminusXC;
        double diffY = direction * edgeTangent(1) / YAminusYC;

        if (abs(XAminusXC) > 1.0e-12 && abs(YAminusYC) > 1.0e-12)
        {
            VectorXd vettX = VectorXd::Ones(polynomialDegree + 2); // ((xa - xc)/diam)^{alphax}
            VectorXd vettY = VectorXd::Ones(polynomialDegree + 2); // ((ya - yc)/diam)^{alphay}
            VectorXd vettDiffX = VectorXd::Ones(polynomialDegree + 2);
            VectorXd vettDiffY = VectorXd::Ones(polynomialDegree + 2);

            for (unsigned int i = 0; i < (polynomialDegree + 1); i++)
            {
                vettX[i + 1] = vettX[i] * X;
                vettY[i + 1] = vettY[i] * Y;
                vettDiffX[i + 1] = vettDiffX[i] * diffX;
                vettDiffY[i + 1] = vettDiffY[i] * diffY;
            }

            Cmatrixkp1[e](0, 0) = 1.0;
            unsigned int offsetRow = 1;

            for (unsigned int N = 1; N <= (polynomialDegree + 1); N++)
            {
                for (unsigned int alphay = 0; alphay <= N; alphay++)
                {
                    for (unsigned int j = 0; j <= alphay; j++)
                    {
                        for (unsigned int i = 0; i <= (N - alphay); i++)
                        {
                            Cmatrixkp1[e](offsetRow, (j + i)) =
                                Cmatrixkp1[e](offsetRow, (j + i)) + BinomialCoefficients[N - alphay][i] * vettDiffX[i] *
                                                                        BinomialCoefficients[alphay][j] * vettDiffY[j];
                        }
                    }
                    Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) =
                        Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) * vettX[N - alphay] * vettY[alphay];
                    offsetRow++;
                }
            }
        }
        else if (abs(XAminusXC) > 1.0e-12 && abs(YAminusYC) <= 1.0e-12)
        {
            VectorXd vettX = VectorXd::Ones(polynomialDegree + 2); // ((xa - xc)/diam)^{alphax}
            VectorXd vettY = VectorXd::Ones(polynomialDegree + 2); // ((yB - yA)/diam)^{alphay}
            VectorXd vettDiffX = VectorXd::Ones(polynomialDegree + 2);

            double YBminusYAForInvDiameter = direction * edgeTangent(1) * invDiameter;

            for (unsigned int i = 0; i < (polynomialDegree + 1); i++)
            {
                vettX[i + 1] = vettX[i] * X;
                vettY[i + 1] = vettY[i] * YBminusYAForInvDiameter;
                vettDiffX[i + 1] = vettDiffX[i] * diffX;
            }

            Cmatrixkp1[e](0, 0) = 1.0;
            unsigned int offsetRow = 1;

            for (unsigned int N = 1; N <= (polynomialDegree + 1); N++)
            {
                for (unsigned int alphay = 0; alphay <= N; alphay++)
                {
                    unsigned int j = alphay;
                    for (unsigned int i = 0; i <= (N - alphay); i++)
                    {
                        Cmatrixkp1[e](offsetRow, (j + i)) =
                            Cmatrixkp1[e](offsetRow, (j + i)) + BinomialCoefficients[N - alphay][i] * vettDiffX[i];
                    }

                    Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) =
                        Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) * vettX[N - alphay] * vettY[alphay];
                    offsetRow++;
                }
            }
        }
        else if (abs(XAminusXC) <= 1.0e-12 && abs(YAminusYC) > 1.0e-12)
        {
            VectorXd vettX = VectorXd::Ones(polynomialDegree + 2); // ((xB - xA)/diam)^{alphax}
            VectorXd vettY = VectorXd::Ones(polynomialDegree + 2); // ((ya - yc)/diam)^{alphay}
            VectorXd vettDiffY = VectorXd::Ones(polynomialDegree + 2);

            double XBminusXAForInvDiameter = direction * edgeTangent(0) * invDiameter;

            for (unsigned int i = 0; i < (polynomialDegree + 1); i++)
            {
                vettX[i + 1] = vettX[i] * XBminusXAForInvDiameter;
                vettY[i + 1] = vettY[i] * Y;
                vettDiffY[i + 1] = vettDiffY[i] * diffY;
            }

            Cmatrixkp1[e](0, 0) = 1.0;
            unsigned int offsetRow = 1;

            for (unsigned int N = 1; N <= (polynomialDegree + 1); N++)
            {
                for (unsigned int alphay = 0; alphay <= N; alphay++)
                {
                    for (unsigned int j = 0; j <= alphay; j++)
                    {
                        unsigned int i = N - alphay;
                        Cmatrixkp1[e](offsetRow, (j + i)) =
                            Cmatrixkp1[e](offsetRow, (j + i)) + BinomialCoefficients[alphay][j] * vettDiffY[j];
                    }
                    Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) =
                        Cmatrixkp1[e].block(offsetRow, 0, 1, order_1D) * vettX[N - alphay] * vettY[alphay];
                    offsetRow++;
                }
            }
        }
        else
        {
            throw std::runtime_error("Cmatrix is wrong");
        }
    }
}
//****************************************************************************
} // namespace MCC
} // namespace VEM
} // namespace Polydim
