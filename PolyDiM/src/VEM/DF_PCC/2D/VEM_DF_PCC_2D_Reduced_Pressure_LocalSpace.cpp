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

#include "VEM_DF_PCC_2D_Reduced_Pressure_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
//****************************************************************************
VEM_DF_PCC_2D_Pressure_LocalSpace_Data VEM_DF_PCC_2D_Reduced_Pressure_LocalSpace::CreateLocalSpace(
    const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
    const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_DF_PCC_2D_Pressure_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    InitializeProjectorsComputation(reference_element_data,
                                    polygon.Centroid,
                                    polygon.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Pressure_LocalSpace::InitializeProjectorsComputation(
    const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
    const Eigen::Vector3d &polygonCentroid,
    const double &polygonDiameter,
    const Eigen::MatrixXd &internalQuadraturePoints,
    const Eigen::VectorXd &internalQuadratureWeights,
    VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace) const
{
    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.Nk = (reference_element_data.Order + 1) * (reference_element_data.Order + 2) / 2;
    localSpace.Nkm1 = reference_element_data.Order * (reference_element_data.Order + 1) / 2;

    // Compute Vandermonde matrices.
    localSpace.Centroid = polygonCentroid;
    localSpace.Diameter = polygonDiameter;
    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polygonCentroid, polygonDiameter);

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;
}
//****************************************************************************
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim
