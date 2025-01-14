#include "VEM_DF_PCC_Utilities.hpp"
#include "lagrange_1D.hpp"

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
template <unsigned short dimension>
Eigen::VectorXd VEM_DF_PCC_Utilities<dimension>::ComputeEdgeBasisCoefficients(const unsigned int &order,
                                                                              const Eigen::VectorXd &edgeInternalPoints) const
{
    // Compute basis function coefficients on the generic edge.
    VectorXd interpolation_points_x(order + 1);
    interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
    return Interpolation::Lagrange::Lagrange_1D_cofficients(interpolation_points_x);
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                                              const unsigned int &order,
                                                              const Eigen::VectorXd &edgeBasisCoefficients,
                                                              const Eigen::VectorXd &pointsCurvilinearCoordinates) const
{
    VectorXd interpolation_points_x(order + 1);
    interpolation_points_x << 0.0, 1.0, edgeInternalPoints.transpose();
    return Interpolation::Lagrange::Lagrange_1D_values(interpolation_points_x, edgeBasisCoefficients, pointsCurvilinearCoordinates);
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeDofiDofiStabilizationMatrix(const std::vector<MatrixXd> &projector,
                                                                             const double &coefficient,
                                                                             const std::vector<Eigen::MatrixXd> &dmatrix) const
{
    Eigen::MatrixXd staBmatrix = dmatrix[0] * projector[0];

    for (unsigned int d = 1; d < dimension; d++)
        staBmatrix += dmatrix[d] * projector[d];

    staBmatrix.diagonal().array() -= 1;

    // staBmatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
    staBmatrix = coefficient * staBmatrix.transpose() * staBmatrix;

    return staBmatrix;
}
//****************************************************************************
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim
