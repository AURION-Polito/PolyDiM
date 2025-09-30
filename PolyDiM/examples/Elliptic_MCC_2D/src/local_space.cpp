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

#include "local_space.hpp"
#include <memory>

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{
namespace local_space
{
//***************************************************************************
ReferenceElement_Data CreateReferenceElement(const Program_configuration::MethodTypes &method_type, const unsigned int method_order)
{
    ReferenceElement_Data reference_element_data;
    reference_element_data.Method_Type = method_type;
    reference_element_data.Order = method_order;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        switch (reference_element_data.Method_Type)
        {
        case Program_configuration::MethodTypes::VEM_MCC:
            reference_element_data.VEM_Type = VEM::MCC::VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_MCC_Partial:
            reference_element_data.VEM_Type = VEM::MCC::VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Partial_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_MCC_Ortho:
            reference_element_data.VEM_Type = VEM::MCC::VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
            reference_element_data.VEM_Type = VEM::MCC::VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_EdgeOrtho_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho:
            reference_element_data.VEM_Type = VEM::MCC::VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace;
            break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        reference_element_data.VEM_ReferenceElement_Pressure =
            Polydim::VEM::MCC::create_VEM_MCC_2D_pressure_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data_Pressure =
            reference_element_data.VEM_ReferenceElement_Pressure->Create(method_order);
        reference_element_data.VEM_ReferenceElement_Velocity =
            Polydim::VEM::MCC::create_VEM_MCC_2D_velocity_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data_Velocity =
            reference_element_data.VEM_ReferenceElement_Velocity->Create(method_order);
        reference_element_data.VEM_LocalSpace_Velocity =
            Polydim::VEM::MCC::create_VEM_MCC_2D_velocity_local_space(reference_element_data.VEM_Type);
        reference_element_data.VEM_LocalSpace_Pressure =
            Polydim::VEM::MCC::create_VEM_MCC_2D_pressure_local_space(reference_element_data.VEM_Type);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
}
//***************************************************************************
LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        local_space_data.VEM_Geometry = {config.GeometricTolerance1D(),
                                         config.GeometricTolerance2D(),
                                         mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsCentroids.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsAreas.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsDiameters.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsTriangulations.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeNormals.at(cell2D_index)};

        local_space_data.VEM_LocalSpace_Data_Velocity =
            reference_element_data.VEM_LocalSpace_Velocity->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data_Velocity,
                                                                             local_space_data.VEM_Geometry);

        local_space_data.VEM_LocalSpace_Data_Pressure =
            reference_element_data.VEM_LocalSpace_Pressure->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data_Pressure,
                                                                             local_space_data.VEM_Geometry);

        local_space_data.VEM_LocalSpace_Data_Pressure.Qmatrix = local_space_data.VEM_LocalSpace_Data_Velocity.QmatrixKp1.topLeftCorner(
            local_space_data.VEM_LocalSpace_Data_Velocity.Nk,
            local_space_data.VEM_LocalSpace_Data_Velocity.Nk);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return local_space_data;
}
//***************************************************************************
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                                          const LocalSpace_Data &local_space_data,
                                                          const Polydim::VEM::MCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return reference_element_data.VEM_LocalSpace_Velocity->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data_Velocity,
                                                                                           projectionType);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return reference_element_data.VEM_LocalSpace_Pressure->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data_Pressure);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const ReferenceElement_Data &reference_element_data,
                                                       const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return reference_element_data.VEM_LocalSpace_Velocity->ComputeBasisFunctionsDivergenceValues(
            local_space_data.VEM_LocalSpace_Data_Velocity);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                                    const ReferenceElement_Data &reference_element_data,
                                                    const LocalSpace_Data &local_space_data,
                                                    const Eigen::MatrixXd &edge_quadrature_points)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho: {
        const double direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        return direction * Eigen::MatrixXd::Identity(reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs1D,
                                                     reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs1D);
    }
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = local_space_data.VEM_Geometry.Tolerance1D;
        geometryUtilitiesConfig.Tolerance2D = local_space_data.VEM_Geometry.Tolerance2D;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        const double direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;
        const unsigned int num_vertices = local_space_data.VEM_Geometry.Vertices.cols();
        const Eigen::Vector3d origin = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index]
                                           ? local_space_data.VEM_Geometry.Vertices.col(edge_local_index)
                                           : local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_vertices);
        const Eigen::Vector3d end = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index]
                                        ? local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_vertices)
                                        : local_space_data.VEM_Geometry.Vertices.col(edge_local_index);

        Eigen::MatrixXd point_curvilinear_coordinates = Eigen::MatrixXd::Zero(3, edge_quadrature_points.cols());
        for (unsigned int p = 0; p < edge_quadrature_points.cols(); p++)
            point_curvilinear_coordinates(0, p) =
                geometryUtilities.PointCurvilinearCoordinate(edge_quadrature_points.col(p), origin, end);

        Utilities::Monomials_1D monomials_1D;
        const Eigen::MatrixXd VanderBoundary1D =
            monomials_1D
                .Vander(reference_element_data.VEM_ReferenceElement_Data_Velocity.Monomials_1D,
                        point_curvilinear_coordinates,
                        Eigen::Vector3d::Constant(0.0),
                        1.0)
                .leftCols(reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1) *
            reference_element_data.VEM_ReferenceElement_Data_Velocity.edge_ortho.QmatrixKp1_1D
                .topLeftCorner(reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1,
                               reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1)
                .transpose();

        return (1.0 / local_space_data.VEM_Geometry.EdgesLength(edge_local_index)) * direction * VanderBoundary1D;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Gedim::Quadrature::QuadratureData EdgeQuadrature(const ReferenceElement_Data &reference_element_data,
                                                 const LocalSpace_Data &local_space_data,
                                                 const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        Gedim::Quadrature::QuadratureData quadrature;
        unsigned int num_quadrature_points =
            reference_element_data.VEM_ReferenceElement_Data_Velocity.Quadrature.ReferenceSegmentQuadrature.Points.cols();
        quadrature.Points = local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Points.middleCols(
            num_quadrature_points * edge_local_index,
            num_quadrature_points);
        quadrature.Weights = local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Weights.segment(
            num_quadrature_points * edge_local_index,
            num_quadrature_points);
        return quadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return local_space_data.VEM_LocalSpace_Data_Velocity.InternalQuadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
unsigned int VelocitySize(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return local_space_data.VEM_LocalSpace_Data_Velocity.NumBasisFunctions;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd StabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const VEM::MCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        return reference_element_data.VEM_LocalSpace_Velocity->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data_Velocity,
                                                                                                  projectionType);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Gedim::Quadrature::QuadratureData EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                                      const LocalSpace_Data &local_space_data,
                                                      const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho: {
        unsigned int num_edge_dofs = reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs1D;

        Gedim::Quadrature::QuadratureData edge_dofs_coordinates;
        edge_dofs_coordinates.Points =
            local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Points.middleCols(num_edge_dofs * edge_local_index,
                                                                                                          num_edge_dofs);
        return edge_dofs_coordinates;
    }
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        Gedim::Quadrature::QuadratureData quadrature;
        unsigned int num_quadrature_points =
            reference_element_data.VEM_ReferenceElement_Data_Velocity.Quadrature.ReferenceSegmentQuadrature.Points.cols();
        quadrature.Points = local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Points.middleCols(
            num_quadrature_points * edge_local_index,
            num_quadrature_points);
        quadrature.Weights = local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Weights.segment(
            num_quadrature_points * edge_local_index,
            num_quadrature_points);
        return quadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::VectorXd EdgeDofs(const ReferenceElement_Data &reference_element_data,
                         const LocalSpace_Data &local_space_data,
                         const unsigned int edge_local_index,
                         const Gedim::Quadrature::QuadratureData &edge_dofs_coordinates,
                         const Eigen::VectorXd &strong_values)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho: {
        const double direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        return direction * strong_values;
    }
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = local_space_data.VEM_Geometry.Tolerance1D;
        geometryUtilitiesConfig.Tolerance2D = local_space_data.VEM_Geometry.Tolerance2D;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        const double direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;
        const unsigned int num_vertices = local_space_data.VEM_Geometry.Vertices.cols();
        const Eigen::Vector3d origin = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index]
                                           ? local_space_data.VEM_Geometry.Vertices.col(edge_local_index)
                                           : local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_vertices);
        const Eigen::Vector3d end = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index]
                                        ? local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_vertices)
                                        : local_space_data.VEM_Geometry.Vertices.col(edge_local_index);

        Eigen::MatrixXd point_curvilinear_coordinates = Eigen::MatrixXd::Zero(3, edge_dofs_coordinates.Points.cols());
        for (unsigned int p = 0; p < edge_dofs_coordinates.Points.cols(); p++)
            point_curvilinear_coordinates(0, p) =
                geometryUtilities.PointCurvilinearCoordinate(edge_dofs_coordinates.Points.col(p), origin, end);

        Utilities::Monomials_1D monomials_1D;
        const Eigen::MatrixXd VanderBoundary1D =
            monomials_1D
                .Vander(reference_element_data.VEM_ReferenceElement_Data_Velocity.Monomials_1D,
                        point_curvilinear_coordinates,
                        Eigen::Vector3d::Constant(0.0),
                        1.0)
                .leftCols(reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1) *
            reference_element_data.VEM_ReferenceElement_Data_Velocity.edge_ortho.QmatrixKp1_1D
                .topLeftCorner(reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1,
                               reference_element_data.VEM_ReferenceElement_Data_Velocity.Order + 1)
                .transpose();

        return direction * VanderBoundary1D.transpose() * edge_dofs_coordinates.Weights.asDiagonal() * strong_values;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::array<std::array<unsigned int, 4>, 2> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        std::array<std::array<unsigned int, 4>, 2> result;
        result[0] = {reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs0D,
                     reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs1D,
                     reference_element_data.VEM_ReferenceElement_Data_Velocity.NumDofs2D,
                     0};
        result[1] = {reference_element_data.VEM_ReferenceElement_Data_Pressure.NumDofs0D,
                     reference_element_data.VEM_ReferenceElement_Data_Pressure.NumDofs1D,
                     reference_element_data.VEM_ReferenceElement_Data_Pressure.NumDofs2D,
                     0};
        return result;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    Performance_Data performance;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::VEM_MCC:
    case Program_configuration::MethodTypes::VEM_MCC_Partial:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho:
    case Program_configuration::MethodTypes::VEM_MCC_EdgeOrtho:
    case Program_configuration::MethodTypes::VEM_MCC_Ortho_EdgeOrtho: {
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        performance.VEM_Performance_Data.Analysis =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(),
                                        reference_element_data.VEM_ReferenceElement_Data_Velocity.MonomialsKp1,
                                        *reference_element_data.VEM_LocalSpace_Velocity,
                                        local_space_data.VEM_LocalSpace_Data_Velocity);

        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data_Velocity.InternalQuadrature.Weights.size();
        performance.VEM_Performance_Data.NumBoundaryQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data_Velocity.BoundaryQuadrature.Quadrature.Weights.size();
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return performance;
}
//***************************************************************************
} // namespace local_space
} // namespace Elliptic_MCC_2D
} // namespace examples
} // namespace Polydim
