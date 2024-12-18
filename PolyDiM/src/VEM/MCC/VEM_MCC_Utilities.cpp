#include "VEM_MCC_Utilities.hpp"

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
MatrixXd VEM_MCC_Utilities<dimension>::ComputeStabilizationMatrix(const MatrixXd &pi0k,
                                                                  const double &measure,
                                                                  const Eigen::MatrixXd &DMatrix) const
{
    MatrixXd stabMatrixPi0k = DMatrix * pi0k;
    stabMatrixPi0k.diagonal().array() -= 1;
    // stabMatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
    stabMatrixPi0k = measure * stabMatrixPi0k.transpose() * stabMatrixPi0k;

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
} // namespace MCC
} // namespace VEM
} // namespace Polydim
