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

#include "FEM_Triangle_RT_MCC_2D_LocalSpace.hpp"
#include "CommonUtilities.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace MCC
{
// ***************************************************************************
FEM_Triangle_RT_MCC_2D_LocalSpace_Data FEM_Triangle_RT_MCC_2D_LocalSpace::CreateLocalSpace(
    const FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
    const FEM_MCC_2D_Polygon_Geometry &polygon) const
{

    if (polygon.Vertices.cols() != 3)
        throw std::runtime_error("The element must be a triangle");

    FEM_Triangle_RT_MCC_2D_LocalSpace_Data localSpace;

    Gedim::MapTriangle mapTriangle;
    localSpace.MapData = mapTriangle.Compute(polygon.Vertices);

    localSpace.Order = reference_element_data.Order;
    localSpace.NumVelocityBasisFunctions = reference_element_data.reference_element_data_velocity.NumBasisFunctions;
    localSpace.NumPressureBasisFunctions = reference_element_data.reference_element_data_pressure.NumBasisFunctions;

    for (unsigned int e = 0; e < 3; e++)
        localSpace.EdgesDirection[e] = polygon.EdgesDirection[e];

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature, localSpace);
    localSpace.BoundaryQuadrature =
        BoundaryQuadrature(reference_element_data, localSpace, reference_element_data.Quadrature.ReferenceSegmentQuadrature, polygon);

    return localSpace;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Triangle_RT_MCC_2D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_Triangle_RT_MCC_2D_LocalSpace_Data &localSpace) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    Gedim::MapTriangle mapTriangle;
    quadrature.Points = Gedim::MapTriangle::F(localSpace.MapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapTriangle::DetJ(localSpace.MapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
std::vector<Gedim::Quadrature::QuadratureData> FEM_Triangle_RT_MCC_2D_LocalSpace::BoundaryQuadrature(
    const FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
    const FEM_Triangle_RT_MCC_2D_LocalSpace_Data &localSpace,
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_MCC_2D_Polygon_Geometry &polygon) const
{
    const unsigned int num_edges = polygon.EdgesDirection.size();
    std::vector<Gedim::Quadrature::QuadratureData> edges_quadrature(num_edges);

    Gedim::MapTriangle mapTriangle;

    for (unsigned int e = 0; e < num_edges; ++e)
    {
        auto &edge_quadrature = edges_quadrature.at(e);
        const double edge_length = polygon.EdgesLength[e];

        edge_quadrature.Points = Gedim::MapTriangle::F(
            localSpace.MapData,
            reference_element_data.BoundaryQuadrature.at(localSpace.EdgesDirection)
                .Quadrature.Points.middleCols(e * reference_element_data.reference_element_data_velocity.NumDofs1D,
                                              reference_element_data.reference_element_data_velocity.NumDofs1D));

        edge_quadrature.Weights = reference_quadrature.Weights * std::abs(edge_length);
    }

    return edges_quadrature;
}
// ***************************************************************************
std::vector<Eigen::MatrixXd> FEM_Triangle_RT_MCC_2D_LocalSpace::MapVelocityValues(
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
    const std::vector<Eigen::MatrixXd> &referenceValues) const
{
    std::vector<Eigen::MatrixXd> velocity_values(2, Eigen::MatrixXd::Zero(referenceValues[0].rows(), referenceValues[0].cols()));
    velocity_values[0] = (1.0 / local_space.MapData.DetB) * (local_space.MapData.B(0, 0) * referenceValues[0] +
                                                             local_space.MapData.B(0, 1) * referenceValues[1]);
    velocity_values[1] = (1.0 / local_space.MapData.DetB) * (local_space.MapData.B(1, 0) * referenceValues[0] +
                                                             local_space.MapData.B(1, 1) * referenceValues[1]);

    return velocity_values;
}
// ***************************************************************************
Eigen::MatrixXd FEM_Triangle_RT_MCC_2D_LocalSpace::MapPressureValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                                                     const Eigen::MatrixXd &referenceValues) const
{
    Gedim::Utilities::Unused(local_space);

    return referenceValues;
}
// ***************************************************************************
Eigen::MatrixXd FEM_Triangle_RT_MCC_2D_LocalSpace::MapVelocityDivergenceValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                                                               const Eigen::MatrixXd &referenceDerivateValues) const
{

    Eigen::MatrixXd divergence_values = (1.0 / local_space.MapData.DetB) * referenceDerivateValues;
    return divergence_values;
}
// ***************************************************************************
} // namespace MCC
} // namespace FEM
} // namespace Polydim
