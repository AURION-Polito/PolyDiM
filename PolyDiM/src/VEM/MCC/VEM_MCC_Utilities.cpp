#include "VEM_MCC_Utilities.hpp"
#include "VEM_MCC_2D_VelocityLocalSpace_Data.hpp"

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
template<unsigned short dimension>
MatrixXd VEM_MCC_Utilities<dimension>::ComputeStabilizationMatrixPi0k(const MatrixXd& pi0k,
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
void VEM_MCC_Utilities<dimension>::ComputeL2Projectors(const unsigned int& numberProjectorBasisFunctions,
                                                       const unsigned int& numberBasisFunctions,
                                                       const unsigned int& numberInternalBasisFunctions,
                                                       const unsigned int& order,
                                                       const MatrixXd& piNabla,
                                                       const MatrixXd& Hmatrix,
                                                       const double& measure,
                                                       const unsigned int& Nkm1,
                                                       const Eigen::LLT<Eigen::MatrixXd>& H_km1_LLT,
                                                       Eigen::MatrixXd& pi0km1,
                                                       Eigen::MatrixXd& pi0k) const
{
    MatrixXd Cmatrix(numberProjectorBasisFunctions, numberBasisFunctions);
    // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
    Cmatrix.bottomRows(numberProjectorBasisFunctions - numberInternalBasisFunctions) =
        Hmatrix.bottomRows(numberProjectorBasisFunctions - numberInternalBasisFunctions)*piNabla;
    if (order > 1)
    {
        Cmatrix.topLeftCorner(numberInternalBasisFunctions,
                              numberBasisFunctions - numberInternalBasisFunctions).setZero();
        // \int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
        Cmatrix.topRightCorner(numberInternalBasisFunctions, numberInternalBasisFunctions) =
            measure*MatrixXd::Identity(numberInternalBasisFunctions, numberInternalBasisFunctions);
    }

    pi0km1 = H_km1_LLT.solve(Cmatrix.topRows(Nkm1));
    pi0k = Hmatrix.llt().solve(Cmatrix);
}
//****************************************************************************
template<unsigned short dimension>
MatrixXd VEM_MCC_Utilities<dimension>::ComputePolynomialBasisDofs(const double& polytopeMeasure,
                                                                  const VEM_MCC_VelocityLocalSpace_Data& localSpace) const
{
    MatrixXd polynomialBasisDofs = MatrixXd::Zero(localSpace.NumBasisFunctions,
                                                  localSpace.Dimension * localSpace.Nk);

    polynomialBasisDofs.topRows(localSpace.NumBoundaryBasisFunctions) = localSpace.GkVanderBoundaryTimesNormal.transpose();

    if(localSpace.Order  > 0)
    {
        polynomialBasisDofs.block(localSpace.NumBoundaryBasisFunctions, 0,
                                  localSpace.NumNablaInternalBasisFunctions, localSpace.Dimension * localSpace.Nk)
            = (1.0 / polygonMeasure) * localSpace.Gmatrix.topRows(localSpace.NumNablaInternalBasisFunctions);

        polynomialBasisDofs.bottomRows(localSpace.NumBigOPlusInternalBasisFunctions)
            = (1.0 / polygonMeasure) * localSpace.Gmatrix.bottomRows(localSpace.NumBigOPlusInternalBasisFunctions);

    }
    return polynomialBasisDofs;
}
//****************************************************************************
}
}
}
