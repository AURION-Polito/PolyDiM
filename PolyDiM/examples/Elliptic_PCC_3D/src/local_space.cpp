#include "local_space.hpp"
#include <memory>

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
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
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        reference_element_data.FEM_ReferenceElement_3D = std::make_unique<FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement>();
        reference_element_data.FEM_ReferenceElement_Data_3D = reference_element_data.FEM_ReferenceElement_3D->Create(method_order);
        reference_element_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_Tetrahedron_PCC_3D_LocalSpace>();
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        switch (reference_element_data.Method_Type)
        {
        case Program_configuration::MethodTypes::VEM_PCC:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_Inertia_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_Ortho_LocalSpace;
            break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        reference_element_data.VEM_ReferenceElement_2D =
            Polydim::VEM::PCC::create_VEM_PCC_3D_reference_element_2D(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data_2D = reference_element_data.VEM_ReferenceElement_2D->Create(method_order);

        reference_element_data.VEM_ReferenceElement_3D =
            Polydim::VEM::PCC::create_VEM_PCC_3D_reference_element_3D(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data_3D = reference_element_data.VEM_ReferenceElement_3D->Create(method_order);

        reference_element_data.VEM_LocalSpace =
            Polydim::VEM::PCC::create_VEM_PCC_3D_local_space_3D(reference_element_data.VEM_Type);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
}
//***************************************************************************
LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                 const unsigned int cell3D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        local_space_data.FEM_Geometry = {config.GeometricTolerance1D(),
                                         config.GeometricTolerance2D(),
                                         config.GeometricTolerance3D(),
                                         mesh_geometric_data.Cell3DsVertices.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsEdges.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFaces.at(cell3D_index),
                                         {},
                                         mesh_geometric_data.Cell3DsEdgeDirections.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesRotationMatrices.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesTranslations.at(cell3D_index)};

        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(cell3D_index).size();
        local_space_data.FEM_Geometry.Faces_2D_Geometry.resize(numFaces);
        for (unsigned int f = 0; f < numFaces; f++)
        {
            local_space_data.FEM_Geometry.Faces_2D_Geometry[f] = {
                mesh_geometric_data.Cell3DsFaces2DVertices.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdgeDirections.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdgeLengths.at(cell3D_index)[f]};
        }

        local_space_data.FEM_LocalSpace_Data =
            reference_element_data.FEM_LocalSpace->CreateLocalSpace(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                    local_space_data.FEM_Geometry);
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(cell3D_index).size();
        for (unsigned int f = 0; f < numFaces; f++)
        {
            local_space_data.VEM_Geometry.PolygonalFaces.push_back(
                {config.GeometricTolerance1D(),
                 config.GeometricTolerance2D(),
                 mesh_geometric_data.Cell3DsFaces2DVertices.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFaces2DCentroids.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesAreas.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesDiameters.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFaces2DTriangulations.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesEdgeLengths.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesEdgeDirections.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(cell3D_index)[f],
                 mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(cell3D_index)[f]});
        }

        local_space_data.VEM_Geometry.Polyhedron = {config.GeometricTolerance1D(),
                                                    config.GeometricTolerance2D(),
                                                    config.GeometricTolerance3D(),
                                                    mesh_geometric_data.Cell3DsVertices.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsEdges.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsFaces.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsCentroids.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsVolumes.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsDiameters.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsTetrahedronPoints.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsFacesRotationMatrices.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsFacesTranslations.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsFacesNormals.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsFacesNormalDirections.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsEdgeDirections.at(cell3D_index),
                                                    mesh_geometric_data.Cell3DsEdgeTangents.at(cell3D_index)};

        local_space_data.VEM_LocalSpace_Data =
            reference_element_data.VEM_LocalSpace->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data_2D,
                                                                    reference_element_data.VEM_ReferenceElement_Data_3D,
                                                                    local_space_data.VEM_Geometry.PolygonalFaces,
                                                                    local_space_data.VEM_Geometry.Polyhedron);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return local_space_data;
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                     const LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                  local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data, projectionType);
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
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                            local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_LocalSpace_Data,
                                                                                            projectionType);
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
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
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
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
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
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return Eigen::MatrixXd::Zero(local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                     local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data,
                                                                                         projectionType);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Gedim::Quadrature::QuadratureData FaceQuadrature(const ReferenceElement_Data &reference_element_data,
                                                 const LocalSpace_Data &local_space_data,
                                                 const unsigned int face_local_index,
                                                 unsigned int &quadrature_offset)
{
    Gedim::Quadrature::QuadratureData faceQuadrature;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature.at(face_local_index);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        const unsigned int num_face_quadrature_points =
            local_space_data.VEM_LocalSpace_Data.facesLocalSpace[face_local_index].InternalQuadrature.Weights.size();
        faceQuadrature.Points =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Points.block(0, quadrature_offset, 3, num_face_quadrature_points);
        ;
        faceQuadrature.Weights =
            local_space_data.VEM_LocalSpace_Data.facesLocalSpace[face_local_index].InternalQuadrature.Weights;
        quadrature_offset += num_face_quadrature_points;
        return faceQuadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValuesOnFace(const unsigned int &face_local_index,
                                           const ReferenceElement_Data &reference_element_data,
                                           const LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &quadrature_points)
{
    // basis function of face including vertices, edges and face evaluated in quadrature_points
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValuesOnFace(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                        local_space_data.FEM_LocalSpace_Data,
                                                                                        face_local_index);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.FaceProjectedBasisFunctionsValues[face_local_index];
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FaceDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                                      const LocalSpace_Data &local_space_data,
                                                      const unsigned int face_local_index,
                                                      unsigned int &quadrature_offset)
{
    Gedim::Quadrature::QuadratureData face_dofs_coordinates;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        const auto &dof_coordinates = local_space_data.FEM_LocalSpace_Data.Dofs;

        const unsigned int cell2DStartingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof2DsIndex.at(face_local_index);
        const unsigned int cell2DEndingLocalIdex = local_space_data.FEM_LocalSpace_Data.Dof2DsIndex.at(face_local_index + 1);
        const unsigned int num_face_dofs = cell2DEndingLocalIdex - cell2DStartingLocalIdex;

        face_dofs_coordinates.Points =
            (num_face_dofs == 0) ? Eigen::MatrixXd(0, 0) : dof_coordinates.block(0, cell2DStartingLocalIdex, 3, num_face_dofs);

        return face_dofs_coordinates;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        const unsigned int num_face_quadrature_points =
            local_space_data.VEM_LocalSpace_Data.facesLocalSpace.at(face_local_index).InternalQuadrature.Weights.size();
        face_dofs_coordinates.Points =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Points.block(0, quadrature_offset, 3, num_face_quadrature_points);

        face_dofs_coordinates.Weights =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Weights.segment(quadrature_offset,
                                                                                               num_face_quadrature_points);
        quadrature_offset += num_face_quadrature_points;

        return face_dofs_coordinates;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
// ***************************************************************************
Eigen::VectorXd FaceDofs(const ReferenceElement_Data &reference_element_data,
                         const LocalSpace_Data &local_space_data,
                         const unsigned int face_local_index,
                         const Eigen::VectorXd &strong_values,
                         const Gedim::Quadrature::QuadratureData &face_dofs_coordinates)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return strong_values;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.FaceScaledMomentsBasis[face_local_index].transpose() *
               face_dofs_coordinates.Weights.asDiagonal() * strong_values;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
// ***************************************************************************
Eigen::MatrixXd EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
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
        const auto &referenceEdgeDOFsPoint = reference_element_data.VEM_ReferenceElement_Data_2D.Quadrature.ReferenceEdgeDOFsInternalPoints;
        const unsigned int num_edge_dofs = referenceEdgeDOFsPoint.cols();

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const Eigen::Vector3d edge_origin = local_space_data.VEM_Geometry.Polyhedron.EdgesDirection.at(edge_local_index)
                                                ? local_space_data.VEM_Geometry.Polyhedron.Vertices.col(
                                                      local_space_data.VEM_Geometry.Polyhedron.Edges(0, edge_local_index))
                                                : local_space_data.VEM_Geometry.Polyhedron.Vertices.col(
                                                      local_space_data.VEM_Geometry.Polyhedron.Edges(1, edge_local_index));

        const Eigen::Vector3d edge_tangent = local_space_data.VEM_Geometry.Polyhedron.EdgesTangent.col(edge_local_index);
        const double edge_direction = local_space_data.VEM_Geometry.Polyhedron.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        Eigen::MatrixXd edge_dofs_coordinates = Eigen::MatrixXd::Zero(3, num_edge_dofs);
        for (unsigned int r = 0; r < num_edge_dofs; r++)
        {
            edge_dofs_coordinates.col(r) << edge_origin + edge_direction * referenceEdgeDOFsPoint(0, r) * edge_tangent;
        }

        return edge_dofs_coordinates;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + "not supported");
    }
}
//***************************************************************************
std::array<unsigned int, 4> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data)
{
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        return {reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs0D,
                reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs1D,
                reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs2D,
                reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs3D};
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return {reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs0D,
                reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs1D,
                reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs2D,
                reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs3D};
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
    case Program_configuration::MethodTypes::FEM_Tetrahedron_PCC: {
        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.FEM_LocalSpace_Data.InternalQuadrature.Weights.size();
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.VEM_Performance_Data.Analysis =
            performanceAnalysis.Compute(Polydim::VEM::Utilities::VEM_Monomials_3D(),
                                        reference_element_data.VEM_ReferenceElement_Data_3D.Monomials,
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
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim
