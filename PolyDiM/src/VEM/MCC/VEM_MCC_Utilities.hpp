#ifndef __VEM_MCC_Utilities_HPP
#define __VEM_MCC_Utilities_HPP

#include "Eigen/Eigen"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
enum struct ProjectionTypes
{
    Pi0k = 1,
};

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

    Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Eigen::MatrixXd &projector,
                                                       const double &coefficient,
                                                       const Eigen::MatrixXd &DMatrix) const;

    void MonomialTraceOnEdges(const unsigned int &polynomialDegree,
                              const Eigen::MatrixXd &polygonVertices,
                              const double &polygonDiameter,
                              const Eigen::Vector3d &polygonCentroid,
                              const std::vector<bool> &edgeDirections,
                              const Eigen::MatrixXd &edgeTangents,
                              std::vector<Eigen::MatrixXd> &Cmatrixkp1) const;

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const Eigen::MatrixXd &projector, const Eigen::MatrixXd &GVander) const;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
