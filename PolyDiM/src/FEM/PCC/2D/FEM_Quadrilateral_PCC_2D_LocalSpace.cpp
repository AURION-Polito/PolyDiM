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

#include "FEM_Quadrilateral_PCC_2D_LocalSpace.hpp"
#include "CommonUtilities.hpp"
#include "GeometryUtilities.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_Quadrilateral_PCC_2D_LocalSpace_Data FEM_Quadrilateral_PCC_2D_LocalSpace::CreateLocalSpace(
    const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
    const FEM_PCC_2D_Polygon_Geometry &polygon) const
{

    if (polygon.Vertices.cols() != 4)
        throw std::runtime_error("The element must be a Quadrilateral");

    FEM_Quadrilateral_PCC_2D_LocalSpace_Data localSpace;

    const double a4 = polygon.Vertices(0, 0) - polygon.Vertices(0, 1) + polygon.Vertices(0, 2) - polygon.Vertices(0, 3);
    const double b4 = polygon.Vertices(1, 0) - polygon.Vertices(1, 1) + polygon.Vertices(1, 2) - polygon.Vertices(1, 3);

    Gedim::GeometryUtilitiesConfig config;
    config.Tolerance1D = polygon.Tolerance1D;
    config.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometry_utilities(config);

    localSpace.Vertices = polygon.Vertices;
    if (geometry_utilities.IsValueZero(a4, polygon.Tolerance1D) && geometry_utilities.IsValueZero(b4, polygon.Tolerance1D))
    {
        Gedim::MapParallelogram MapParallelogram;
        localSpace.MapData = MapParallelogram.Compute(polygon.Vertices);
        localSpace.B_lap = localSpace.MapData.BInv * localSpace.MapData.BInv.transpose();
        localSpace.quadrilateral_type = QuadrilateralType::Parallelogram;
    }
    else
        localSpace.quadrilateral_type = QuadrilateralType::Generic;

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);

    unsigned int dofCounter = 0;
    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 4; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[dofCounter] = d;
            dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[4]);
    for (unsigned int e = 0; e < 4; e++)
    {
        localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

        const unsigned int edge_origin_index = e;
        const unsigned int edge_end_index = (e + 1) % 4;

        const auto &reference_edge = reference_element_data.Edges_by_vertices.at({edge_origin_index, edge_end_index});

        if (polygon.EdgesDirection.at(e) == reference_edge.second)
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e]; d < localSpace.Dof1DsIndex[e + 1]; d++)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e + 1] - 1; d < UINT_MAX && d >= localSpace.Dof1DsIndex[e]; d--)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
    }

    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[4]);
    localSpace.Dof2DsIndex[1] = localSpace.Dof2DsIndex[0] + reference_element_data.NumDofs2D;

    for (unsigned int d = localSpace.Dof2DsIndex[0]; d < localSpace.Dof2DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

    localSpace.dofs_permutation.resize(localSpace.NumberOfBasisFunctions);
    for (unsigned int d = 0; d < localSpace.NumberOfBasisFunctions; ++d)
        localSpace.dofs_permutation.indices()[localSpace.DofsMeshOrder.at(d)] = d;

    // reorder basis function values with mesh order
    switch (localSpace.quadrilateral_type)
    {
    case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
        localSpace.Dofs =
            MapValues(localSpace, Gedim::MapParallelogram::F(localSpace.MapData, reference_element_data.DofPositions));
    }
    break;
    case Polydim::FEM::PCC::QuadrilateralType::Generic: {
        Gedim::MapQuadrilateral mapQuadrilateral;
        localSpace.Dofs = MapValues(localSpace, mapQuadrilateral.F(localSpace.Vertices, reference_element_data.DofPositions));
    }
    break;
    default:
        throw std::runtime_error("not valid quadrilateral");
    }

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceSquareQuadrature, localSpace);
    localSpace.BoundaryQuadrature =
        BoundaryQuadrature(reference_element_data.BoundaryReferenceElement_Data.ReferenceSegmentQuadrature, polygon);

    return localSpace;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_Quadrilateral_PCC_2D_LocalSpace::MapDerivativeValues(const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                                                               const std::vector<Eigen::MatrixXd> &referenceDerivateValues,
                                                                               const Eigen::MatrixXd &referencePoints) const
{
    switch (local_space.quadrilateral_type)
    {
    case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
        std::vector<Eigen::MatrixXd> basis_functions_reordered(2);
        basis_functions_reordered.at(0).noalias() = local_space.MapData.BInv(0, 0) * referenceDerivateValues.at(0);
        basis_functions_reordered.at(0).noalias() += local_space.MapData.BInv(1, 0) * referenceDerivateValues.at(1);

        basis_functions_reordered.at(1).noalias() = local_space.MapData.BInv(0, 1) * referenceDerivateValues.at(0);
        basis_functions_reordered.at(1).noalias() += local_space.MapData.BInv(1, 1) * referenceDerivateValues.at(1);

        basis_functions_reordered.at(0) *= local_space.dofs_permutation;
        basis_functions_reordered.at(1) *= local_space.dofs_permutation;

        return basis_functions_reordered;
    }
    case Polydim::FEM::PCC::QuadrilateralType::Generic: {
        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
            2,
            Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

        Gedim::MapQuadrilateral mapQuadrilateral;
        const Eigen::MatrixXd invJac = mapQuadrilateral.JInv(local_space.Vertices, referencePoints);
        for (unsigned int p = 0; p < referencePoints.cols(); p++)
        {
            const Eigen::MatrixXd BInv = invJac.block(0, 3 * p, 2, 2);
            for (unsigned int i = 0; i < 2; i++)
            {
                basisFunctionsDerivativeValues[i].row(p) = BInv(i, i) * referenceDerivateValues[i].row(p);
                for (unsigned int j = 0; j < i; j++)
                {
                    basisFunctionsDerivativeValues[i].row(p) += BInv(j, i) * referenceDerivateValues[j].row(p);
                    basisFunctionsDerivativeValues[j].row(p) += BInv(i, j) * referenceDerivateValues[i].row(p);
                }
            }
        }

        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
            2,
            Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

        for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
        {
            const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
            basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
            basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
        }
        return basisFunctionsDerivativeValuesOrdered;
    }
    default:
        throw std::runtime_error("not valid quadrilateral");
    }
}
// ***************************************************************************
Eigen::MatrixXd FEM_Quadrilateral_PCC_2D_LocalSpace::MapLaplacianValues(const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                                                        const std::array<Eigen::MatrixXd, 4> &referenceSecondDerivateValues,
                                                                        const Eigen::MatrixXd &referencePoints) const
{
    Gedim::Utilities::Unused(referencePoints);

    Eigen::MatrixXd laplacian_values =
        Eigen::MatrixXd::Zero(referenceSecondDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions);

    switch (local_space.quadrilateral_type)
    {
    case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
        laplacian_values = local_space.B_lap(0, 0) * referenceSecondDerivateValues.at(0) +
                           local_space.B_lap(0, 1) * referenceSecondDerivateValues.at(1) +
                           local_space.B_lap(1, 0) * referenceSecondDerivateValues.at(2) +
                           local_space.B_lap(1, 1) * referenceSecondDerivateValues.at(3);
    }
    break;
    default:
        throw std::runtime_error("not valid quadrilateral");
    }

    return laplacian_values;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Quadrilateral_PCC_2D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    switch (local_space.quadrilateral_type)
    {
    case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
        quadrature.Points = Gedim::MapParallelogram::F(local_space.MapData, reference_quadrature.Points);
        quadrature.Weights = reference_quadrature.Weights.array() *
                             Gedim::MapParallelogram::DetJ(local_space.MapData, reference_quadrature.Points).array().abs();
    }
    break;
    case Polydim::FEM::PCC::QuadrilateralType::Generic: {
        Gedim::MapQuadrilateral mapQuadrilateral;
        quadrature.Points = mapQuadrilateral.F(local_space.Vertices, reference_quadrature.Points);
        quadrature.Weights = reference_quadrature.Weights.array() *
                             mapQuadrilateral.DetJ(local_space.Vertices, reference_quadrature.Points).array().abs();
    }
    break;
    default:
        throw std::runtime_error("not valid quadrilateral");
    }

    return quadrature;
}
// ***************************************************************************
std::vector<Gedim::Quadrature::QuadratureData> FEM_Quadrilateral_PCC_2D_LocalSpace::BoundaryQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_PCC_2D_Polygon_Geometry &polygon) const
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
} // namespace PCC

} // namespace FEM
} // namespace Polydim
