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
namespace Elastic_PCC_2D
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
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        reference_element_data.FEM_ReferenceElement = std::make_unique<FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement>();
        reference_element_data.FEM_ReferenceElement_Data = reference_element_data.FEM_ReferenceElement->Create(method_order);
        reference_element_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace>();
        reference_element_data.Dimension = reference_element_data.FEM_ReferenceElement_Data.Dimension;
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        switch (reference_element_data.Method_Type)
        {
        case Program_configuration::MethodTypes::VEM_PCC:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Ortho_LocalSpace;
            break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        reference_element_data.VEM_ReferenceElement =
            Polydim::VEM::PCC::create_VEM_PCC_2D_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data = reference_element_data.VEM_ReferenceElement->Create(method_order);
        reference_element_data.VEM_LocalSpace = Polydim::VEM::PCC::create_VEM_PCC_2D_local_space(reference_element_data.VEM_Type);
        reference_element_data.Dimension = reference_element_data.VEM_ReferenceElement_Data.Dimension;
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
}
//***************************************************************************
LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        local_space_data.FEM_Geometry = {config.GeometricTolerance1D(),
                                         config.GeometricTolerance2D(),
                                         mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index)};

        local_space_data.FEM_LocalSpace_Data =
            reference_element_data.FEM_LocalSpace->CreateLocalSpace(reference_element_data.FEM_ReferenceElement_Data,
                                                                    local_space_data.FEM_Geometry);
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
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
        local_space_data.VEM_LocalSpace_Data =
            reference_element_data.VEM_LocalSpace->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data,
                                                                    local_space_data.VEM_Geometry);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return local_space_data;
}
//***************************************************************************
std::vector<Eigen::MatrixXd> BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                                  const LocalSpace_Data &local_space_data,
                                                  const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{

    Eigen::MatrixXd basis_functions_values;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        basis_functions_values =
            reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                               local_space_data.FEM_LocalSpace_Data);
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        basis_functions_values =
            reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data, projectionType);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    const unsigned int num_local_dofs_component = basis_functions_values.cols();
    std::vector<Eigen::MatrixXd> basis_functions_displacement_values(
        reference_element_data.Dimension,
        Eigen::MatrixXd::Zero(basis_functions_values.rows(), reference_element_data.Dimension * num_local_dofs_component));
    for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
        basis_functions_displacement_values[d].middleCols(d * num_local_dofs_component, num_local_dofs_component) = basis_functions_values;

    return basis_functions_displacement_values;
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsLaplacianValues(const ReferenceElement_Data &reference_element_data,
                                              const LocalSpace_Data &local_space_data,
                                              const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsLaplacianValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                           local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsLaplacianValues(local_space_data.VEM_LocalSpace_Data,
                                                                                           projectionType);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                            const ReferenceElement_Data &reference_element_data,
                                            const LocalSpace_Data &local_space_data,
                                            const Eigen::MatrixXd &pointsCurvilinearCoordinates)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        VEM::PCC::VEM_PCC_Utilities<2> utilities;

        const auto &dof_coordinates = local_space_data.FEM_LocalSpace_Data.Dofs;

        const unsigned int cell1DStartingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof1DsIndex.at(edge_local_index);
        const unsigned int cell1DEndingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof1DsIndex.at(edge_local_index + 1);
        const unsigned int num_edge_dofs = cell1DEndingLocalIdex - cell1DStartingLocalIdex;

        const Eigen::VectorXd edgeInternalPoints = Eigen::VectorXd::LinSpaced(num_edge_dofs + 2, 0.0, 1.0).segment(1, num_edge_dofs);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(), reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeValuesOnEdge(reference_element_data.VEM_ReferenceElement_Data,
                                                                          pointsCurvilinearCoordinates);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(const ReferenceElement_Data &reference_element_data,
                                                            const LocalSpace_Data &local_space_data,
                                                            const VEM::PCC::ProjectionTypes &projectionType)
{

    std::vector<Eigen::MatrixXd> basis_functions_derivatives_values;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        basis_functions_derivatives_values = reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(
            reference_element_data.FEM_ReferenceElement_Data,
            local_space_data.FEM_LocalSpace_Data);
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        basis_functions_derivatives_values =
            reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_LocalSpace_Data,
                                                                                         projectionType);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
    const unsigned int num_local_dofs_component = basis_functions_derivatives_values.at(0).cols();
    std::vector<Eigen::MatrixXd> basis_functions_derivatives_displacement_values(
        reference_element_data.Dimension * reference_element_data.Dimension,
        Eigen::MatrixXd::Zero(basis_functions_derivatives_values.at(0).rows(),
                              reference_element_data.Dimension * num_local_dofs_component));
    for (unsigned int d1 = 0; d1 < reference_element_data.Dimension; d1++)
    {
        for (unsigned int d2 = 0; d2 < reference_element_data.Dimension; d2++)
            basis_functions_derivatives_displacement_values[2 * d1 + d2].middleCols(d1 * num_local_dofs_component,
                                                                                    num_local_dofs_component) =
                basis_functions_derivatives_values[d2];
    }

    return basis_functions_derivatives_displacement_values;
}
//***************************************************************************
Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        return local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.InternalQuadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
unsigned int Size(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        return local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.NumBasisFunctions;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd StabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        return Eigen::MatrixXd::Zero(reference_element_data.Dimension * local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                     reference_element_data.Dimension * local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        const Eigen::MatrixXd scalar_stab_matrix =
            reference_element_data.VEM_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data,
                                                                                      projectionType);

        const unsigned int num_local_dofs_component = scalar_stab_matrix.rows();
        Eigen::MatrixXd stab_matrix = Eigen::MatrixXd::Zero(reference_element_data.Dimension * num_local_dofs_component,
                                                            reference_element_data.Dimension * num_local_dofs_component);

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
            stab_matrix.block(d * num_local_dofs_component, d * num_local_dofs_component, num_local_dofs_component, num_local_dofs_component) =
                scalar_stab_matrix;

        return stab_matrix;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        const auto &dof_coordinates = local_space_data.FEM_LocalSpace_Data.Dofs;

        const unsigned int cell1DStartingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof1DsIndex.at(edge_local_index);
        const unsigned int cell1DEndingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof1DsIndex.at(edge_local_index + 1);
        const unsigned int num_edge_dofs = cell1DEndingLocalIdex - cell1DStartingLocalIdex;

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const Eigen::MatrixXd edge_dofs_coordinates = dof_coordinates.block(0, cell1DStartingLocalIdex, 3, num_edge_dofs);

        return edge_dofs_coordinates;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        const auto &referenceEdgeDOFsPoint = reference_element_data.VEM_ReferenceElement_Data.Quadrature.ReferenceEdgeDOFsInternalPoints;
        const unsigned int num_edge_dofs = referenceEdgeDOFsPoint.cols();

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const unsigned int num_edges = local_space_data.VEM_Geometry.Vertices.cols();
        const Eigen::Vector3d edge_origin = local_space_data.VEM_Geometry.EdgesDirection.at(edge_local_index)
                                                ? local_space_data.VEM_Geometry.Vertices.col(edge_local_index)
                                                : local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_edges);

        const Eigen::Vector3d edge_tangent = local_space_data.VEM_Geometry.EdgesTangent.col(edge_local_index);
        const double edge_direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        Eigen::MatrixXd edge_dofs_coordinates = Eigen::MatrixXd::Zero(3, num_edge_dofs);
        for (unsigned int r = 0; r < num_edge_dofs; r++)
        {
            edge_dofs_coordinates.col(r) << edge_origin + edge_direction * referenceEdgeDOFsPoint(0, r) * edge_tangent;
        }

        return edge_dofs_coordinates;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::array<unsigned int, 4> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        return {reference_element_data.FEM_ReferenceElement_Data.NumDofs0D,
                reference_element_data.FEM_ReferenceElement_Data.NumDofs1D,
                reference_element_data.FEM_ReferenceElement_Data.NumDofs2D,
                0};
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return {reference_element_data.VEM_ReferenceElement_Data.NumDofs0D,
                reference_element_data.VEM_ReferenceElement_Data.NumDofs1D,
                reference_element_data.VEM_ReferenceElement_Data.NumDofs2D,
                0};
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
    case Program_configuration::MethodTypes::FEM_Triangle_PCC: {
        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.FEM_LocalSpace_Data.InternalQuadrature.Weights.size();
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.VEM_Performance_Data.Analysis =
            performanceAnalysis.Compute(Polydim::VEM::Utilities::VEM_Monomials_2D(),
                                        reference_element_data.VEM_ReferenceElement_Data.Monomials,
                                        *reference_element_data.VEM_LocalSpace,
                                        local_space_data.VEM_LocalSpace_Data);

        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.InternalQuadrature.Weights.size();
        performance.VEM_Performance_Data.NumBoundaryQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Weights.size();
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return performance;
}
//***************************************************************************
} // namespace local_space
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim
