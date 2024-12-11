#include "VEM_DF_PCC_Utilities.hpp"

using namespace Eigen;
using namespace std;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
template struct VEM_DF_PCC_Utilities<2>;
template struct VEM_DF_PCC_Utilities<3>;
//****************************************************************************
template<unsigned short dimension>
Eigen::VectorXd VEM_DF_PCC_Utilities<dimension>::ComputeEdgeBasisCoefficients(const unsigned int& order,
                                                                           const Eigen::VectorXd& edgeInternalPoints) const
{
    // Compute basis function coefficients on the generic edge.
    const unsigned int numPoints = order+1;
    // build vector of points on the reference edge.
    VectorXd refEdgePoints(numPoints);
    refEdgePoints << 0.0, 1.0, edgeInternalPoints;
    // vector below is defined as differences[i*(numPoints-1)+j] = x_j - x_i ,
    // where x_k is the k-th quadrature point, for any j != i
    MatrixXd differences(numPoints, numPoints-1);
    for(unsigned int i = 0; i < numPoints; ++i)
    {
        unsigned int col = 0;
        for(unsigned int j = 0; j < numPoints; ++j)
        {
            if( j != i)
            {
                differences(i, col) = refEdgePoints(i) - refEdgePoints(j);
                col++;
            }
        }
    }
    // Compute edgeBasisCoefficients[i]:
    // - compute prod_{j != i} (x_j - x_i)
    VectorXd edgeBasisCoefficients = differences.rowwise().prod();
    for(unsigned int i = 0; i < numPoints; ++i)
    {
        // - invert result.
        edgeBasisCoefficients[i] = 1.0 / edgeBasisCoefficients[i];
    }

    return edgeBasisCoefficients;
}
//****************************************************************************
template<unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeStabilizationMatrix(const MatrixXd& piNabla,
                                                                  const double& diameter,
                                                                  const MatrixXd& Dmatrix) const
{
    MatrixXd stabMatrix = Dmatrix * piNabla;
    stabMatrix.diagonal().array() -= 1;
    // stabMatrix = (\Pi^{\nabla,dofs}_order - I)^T * (\Pi^{\nabla,dofs}_order - I).
    stabMatrix = stabMatrix.transpose() * stabMatrix;

    if (dimension == 3)
        stabMatrix *= diameter;

    return stabMatrix;
}
//****************************************************************************
template<unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeStabilizationMatrixPi0k(const MatrixXd& pi0k,
                                                                      const double& measure,
                                                                      const Eigen::MatrixXd& DMatrix) const
{
    MatrixXd stabMatrixPi0k = DMatrix * pi0k;
    stabMatrixPi0k.diagonal().array() -= 1;
    // stabMatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
    stabMatrixPi0k = measure * stabMatrixPi0k.transpose() * stabMatrixPi0k;

    return stabMatrixPi0k;
}
//****************************************************************************
template<unsigned short dimension>
void VEM_DF_PCC_Utilities<dimension>::ComputeL2Projectors(const double& measure,
                                                       const unsigned int& order,
                                                       const unsigned int& Nkm1,
                                                       const unsigned int& Nk,
                                                       const unsigned int& NumInternalBasisFunctions,
                                                       const unsigned int& NumBasisFunctions,
                                                       const Eigen::MatrixXd& Hmatrix,
                                                       const Eigen::MatrixXd& PiNabla,
                                                       MatrixXd &Cmatrix,
                                                       MatrixXd &Pi0km1,
                                                       MatrixXd &Pi0k) const
{
    Cmatrix  = Eigen::MatrixXd::Zero(Nk, NumBasisFunctions);
    // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
    Cmatrix.bottomRows(Nk - NumInternalBasisFunctions) =
        Hmatrix.bottomRows(Nk - NumInternalBasisFunctions) *
        PiNabla;

    if (order > 1)
    {
        Cmatrix.topLeftCorner(NumInternalBasisFunctions,
                              NumBasisFunctions - NumInternalBasisFunctions).setZero();
        //\int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
        Cmatrix.topRightCorner(NumInternalBasisFunctions, NumInternalBasisFunctions) =
            measure * Eigen::MatrixXd::Identity(NumInternalBasisFunctions,
                                                NumInternalBasisFunctions);
    }

    Pi0km1 = Hmatrix.topLeftCorner(Nkm1, Nkm1).llt().solve(Cmatrix.topRows(Nkm1));
    Pi0k = Hmatrix.llt().solve(Cmatrix);
}
//****************************************************************************
template<unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeValuesOnEdge(const Eigen::RowVectorXd& edgeInternalPoints,
                                                           const unsigned int& order,
                                                           const Eigen::VectorXd& edgeBasisCoefficients,
                                                           const Eigen::VectorXd& pointsCurvilinearCoordinates) const
{
    const unsigned int numPoints = pointsCurvilinearCoordinates.size();
    const unsigned int numBasisFunctions = order + 1;
    Array<double, 1, Eigen::Dynamic> refEdgePoints(numBasisFunctions);
    refEdgePoints << 0.0, 1.0, edgeInternalPoints;
    MatrixXd differences(numPoints,
                         numBasisFunctions);
    for(unsigned int i = 0; i < numPoints; ++i)
        differences.row(i) = pointsCurvilinearCoordinates(i) - refEdgePoints;

    MatrixXd values = MatrixXd::Zero(numPoints, numBasisFunctions);

    for (unsigned int i = 0; i < numBasisFunctions; ++i)
    {
        values.col(i).setConstant(edgeBasisCoefficients[i]);
        for (unsigned int j = 0; j < numBasisFunctions; ++j)
        {
            if (i != j)
                values.col(i) = values.col(i).cwiseProduct(differences.col(j));
        }
    }

    return values;
}
//****************************************************************************
}
}
}
