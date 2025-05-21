// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "VEM_PCC_Utilities.hpp"
#include "lagrange_1D.hpp"

using namespace Eigen;
using namespace std;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
template struct VEM_PCC_Utilities<2>;
template struct VEM_PCC_Utilities<3>;
//****************************************************************************
template <unsigned short dimension>
Eigen::VectorXd VEM_PCC_Utilities<dimension>::ComputeEdgeBasisCoefficients(const unsigned int &order,
                                                                           const Eigen::VectorXd &edgeInternalPoints) const
{
    // Compute basis function coefficients on the generic edge.
    VectorXd interpolation_points_x(order + 1);
    interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
    return Interpolation::Lagrange::Lagrange_1D_coefficients(interpolation_points_x);
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_PCC_Utilities<dimension>::ComputeDofiDofiStabilizationMatrix(const MatrixXd &projector,
                                                                          const double &coefficient,
                                                                          const MatrixXd &Dmatrix) const
{
    MatrixXd stabMatrix = Dmatrix * projector;
    stabMatrix.diagonal().array() -= 1;
    // stabMatrix = (\Pi^{\nabla,dofs}_order - I)^T * (\Pi^{\nabla,dofs}_order - I).
    stabMatrix = coefficient * stabMatrix.transpose() * stabMatrix;

    return stabMatrix;
}
//****************************************************************************
template <unsigned short dimension>
void VEM_PCC_Utilities<dimension>::ComputeL2Projectors(const double &measure,
                                                       const unsigned int &order,
                                                       const unsigned int &Nkm1,
                                                       const unsigned int &Nk,
                                                       const unsigned int &NumInternalBasisFunctions,
                                                       const unsigned int &NumBasisFunctions,
                                                       const Eigen::MatrixXd &Hmatrix,
                                                       const Eigen::MatrixXd &PiNabla,
                                                       MatrixXd &Cmatrix,
                                                       MatrixXd &Pi0km1,
                                                       MatrixXd &Pi0k) const
{
    Cmatrix = Eigen::MatrixXd::Zero(Nk, NumBasisFunctions);
    // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
    Cmatrix.bottomRows(Nk - NumInternalBasisFunctions) = Hmatrix.bottomRows(Nk - NumInternalBasisFunctions) * PiNabla;

    if (order > 1)
    {
        Cmatrix.topLeftCorner(NumInternalBasisFunctions, NumBasisFunctions - NumInternalBasisFunctions).setZero();
        //\int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
        Cmatrix.topRightCorner(NumInternalBasisFunctions, NumInternalBasisFunctions) =
            measure * Eigen::MatrixXd::Identity(NumInternalBasisFunctions, NumInternalBasisFunctions);
    }

    Pi0km1 = Hmatrix.topLeftCorner(Nkm1, Nkm1).llt().solve(Cmatrix.topRows(Nkm1));
    Pi0k = Hmatrix.llt().solve(Cmatrix);
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_PCC_Utilities<dimension>::ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                                           const unsigned int &order,
                                                           const Eigen::VectorXd &edgeBasisCoefficients,
                                                           const Eigen::VectorXd &pointsCurvilinearCoordinates) const
{
    VectorXd interpolation_points_x(order + 1);
    interpolation_points_x << 0.0, 1.0, edgeInternalPoints.transpose();
    return Interpolation::Lagrange::Lagrange_1D_values(interpolation_points_x, edgeBasisCoefficients, pointsCurvilinearCoordinates);
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
