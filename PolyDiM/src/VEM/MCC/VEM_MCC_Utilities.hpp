#ifndef __VEM_MCC_Utilities_HPP
#define __VEM_MCC_Utilities_HPP

#include "Eigen/Eigen"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
/// \brief Base class for computing values of basis functions of Mixed Conforming Constant degree
/// Virtual Element Methods.
/// \copyright See top level LICENSE file for details.
template <unsigned short dimension> struct VEM_MCC_Utilities final
{
    Eigen::MatrixXd ComputePolynomialBasisDofs(const double &polytopeMeasure,
                                               const unsigned int &order,
                                               const unsigned int &Nk,
                                               const unsigned int &NumBoundaryBasisFunctions,
                                               const unsigned int &NumNablaInternalBasisFunctions,
                                               const unsigned int &NumBigOPlusInternalBasisFunctions,
                                               const unsigned int &NumBasisFunctions,
                                               const Eigen::MatrixXd &GkVanderBoundaryTimesNormal,
                                               const Eigen::MatrixXd &Gmatrix) const;

    Eigen::MatrixXd ComputeStabilizationMatrix(const Eigen::MatrixXd &pi0k,
                                               const double &measure,
                                               const Eigen::MatrixXd &DMatrix) const;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
