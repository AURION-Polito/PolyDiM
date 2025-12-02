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
    localSpace.B_lap = localSpace.MapData.BInv * localSpace.MapData.BInv.transpose();

    localSpace.Order = reference_element_data.Order;

    return localSpace;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Triangle_RT_MCC_2D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const Gedim::MapTriangle::MapTriangleData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    Gedim::MapTriangle mapTriangle;
    quadrature.Points = Gedim::MapTriangle::F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapTriangle::DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
std::vector<Gedim::Quadrature::QuadratureData> FEM_Triangle_RT_MCC_2D_LocalSpace::BoundaryQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_MCC_2D_Polygon_Geometry &polygon) const
{
    const unsigned int num_edges = polygon.EdgesDirection.size();
    std::vector<Gedim::Quadrature::QuadratureData> edges_quadrature(num_edges);

    for (unsigned int e = 0; e < num_edges; ++e)
    {
        auto &edge_quadrature = edges_quadrature.at(e);
        const bool edge_direction = polygon.EdgesDirection.at(e);
        const double edge_length = polygon.EdgesLength[e];
        const Eigen::Vector3d edge_origin = edge_direction ? polygon.Vertices.col(e) : polygon.Vertices.col((e + 1) % num_edges);

        const Eigen::Vector3d edge_tangent = edge_direction ? +1.0 * polygon.EdgesTangent.col(e)
                                                            : -1.0 * polygon.EdgesTangent.col(e);

        const unsigned int num_quadrature_points = reference_quadrature.Points.cols();
        edge_quadrature.Points.resize(3, num_quadrature_points);
        for (unsigned int q = 0; q < num_quadrature_points; q++)
            edge_quadrature.Points.col(q) = edge_origin + reference_quadrature.Points(0, q) * edge_tangent;

        edge_quadrature.Weights = reference_quadrature.Weights * std::abs(edge_length);
    }

    return edges_quadrature;
}
// ***************************************************************************
} // namespace MCC

} // namespace FEM
} // namespace Polydim
