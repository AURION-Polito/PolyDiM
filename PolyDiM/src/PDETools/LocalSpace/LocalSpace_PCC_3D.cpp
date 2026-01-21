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

#include "LocalSpace_PCC_3D.hpp"
#include "CommonUtilities.hpp"
#include <memory>

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_PCC_3D
{
//***************************************************************************
ReferenceElement_Data CreateReferenceElement(const MethodTypes &method_type, const unsigned int method_order)
{
    ReferenceElement_Data reference_element_data;
    reference_element_data.Method_Type = method_type;
    reference_element_data.Order = method_order;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        reference_element_data.FEM_ReferenceElement_3D = std::make_unique<FEM::PCC::FEM_PCC_3D_ReferenceElement>();
        reference_element_data.FEM_ReferenceElement_Data_3D = reference_element_data.FEM_ReferenceElement_3D->Create(method_order);
        reference_element_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_PCC_3D_LocalSpace>();
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::VEM_PCC:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_LocalSpace;
            break;
        case MethodTypes::VEM_PCC_Inertia:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_Inertia_LocalSpace;
            break;
        case MethodTypes::VEM_PCC_Ortho:
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
LocalSpace_Data CreateLocalSpace(const double &geometric_tolerance_1D,
                                 const double &geometric_tolerance_2D,
                                 const double &geometric_tolerance_3D,
                                 const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                 const unsigned int cell3D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        local_space_data.FEM_Geometry = {geometric_tolerance_1D,
                                         geometric_tolerance_2D,
                                         geometric_tolerance_3D,
                                         mesh_geometric_data.Cell3DsVertices.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsEdges.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFaces.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsEdgeDirections.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesRotationMatrices.at(cell3D_index),
                                         mesh_geometric_data.Cell3DsFacesTranslations.at(cell3D_index),
                                         {}};

        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(cell3D_index).size();
        local_space_data.FEM_Geometry.Faces_2D_Geometry.resize(numFaces);
        for (unsigned int f = 0; f < numFaces; f++)
        {
            local_space_data.FEM_Geometry.Faces_2D_Geometry[f] = {
                geometric_tolerance_1D,
                geometric_tolerance_2D,
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
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {

        local_space_data.VEM_Geometry = {geometric_tolerance_1D,
                                         geometric_tolerance_2D,
                                         geometric_tolerance_3D,
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
                                         mesh_geometric_data.Cell3DsEdgeTangents.at(cell3D_index),
                                         {}};

        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(cell3D_index).size();
        local_space_data.VEM_Geometry.Faces_2D_Geometry.resize(numFaces);
        for (unsigned int f = 0; f < numFaces; f++)
        {
            local_space_data.VEM_Geometry.Faces_2D_Geometry[f] = {
                geometric_tolerance_1D,
                geometric_tolerance_2D,
                mesh_geometric_data.Cell3DsFaces2DVertices.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFaces2DCentroids.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesAreas.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesDiameters.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFaces2DTriangulations.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdgeLengths.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdgeDirections.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(cell3D_index)[f],
                mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(cell3D_index)[f]};
        }

        local_space_data.VEM_LocalSpace_Data =
            reference_element_data.VEM_LocalSpace->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data_2D,
                                                                    reference_element_data.VEM_ReferenceElement_Data_3D,
                                                                    local_space_data.VEM_Geometry);
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
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                  local_space_data.FEM_LocalSpace_Data);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                            local_space_data.FEM_LocalSpace_Data);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return Eigen::MatrixXd::Zero(local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                     local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature.at(face_local_index);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    Gedim::Utilities::Unused(quadrature_points);

    // basis function of face including vertices, edges and face evaluated in quadrature_points
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValuesOnFace(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                                        local_space_data.FEM_LocalSpace_Data,
                                                                                        face_local_index);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {

        face_dofs_coordinates.Points =
            reference_element_data.FEM_LocalSpace->FaceDOFsCoordinates(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                       local_space_data.FEM_LocalSpace_Data,
                                                                       face_local_index);
        return face_dofs_coordinates;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return strong_values;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
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
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->EdgeDOFsCoordinates(reference_element_data.FEM_ReferenceElement_Data_3D,
                                                                          local_space_data.FEM_LocalSpace_Data,
                                                                          edge_local_index);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->EdgeDOFsCoordinates(local_space_data.VEM_LocalSpace_Data, edge_local_index);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + "not supported");
    }
}
//***************************************************************************
PDETools::DOFs::DOFsManager::MeshDOFsInfo SetMeshDOFsInfo(
    const ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> &boundary_info)
{
    PDETools::DOFs::DOFsManager::MeshDOFsInfo mesh_dof_info;

    const unsigned int numCell0Ds = mesh.Cell0DTotalNumber();

    mesh_dof_info.CellsNumDOFs[0].resize(numCell0Ds);
    mesh_dof_info.CellsBoundaryInfo[0].resize(numCell0Ds);

    for (unsigned int c = 0; c < numCell0Ds; ++c)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs0D;
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs0D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        mesh_dof_info.CellsBoundaryInfo[0][c] = boundary_info.at(mesh.Cell0DMarker(c));
    }

    const unsigned int numCell1Ds = mesh.Cell1DTotalNumber();

    mesh_dof_info.CellsNumDOFs[1].resize(numCell1Ds);
    mesh_dof_info.CellsBoundaryInfo[1].resize(numCell1Ds);

    for (unsigned int c = 0; c < numCell1Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.FEM_ReferenceElement_Data_3D.NumDofs1D;
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs1D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
        mesh_dof_info.CellsBoundaryInfo[1][c] = boundary_info.at(mesh.Cell1DMarker(c));
    }

    const unsigned int numCell2Ds = mesh.Cell2DTotalNumber();

    mesh_dof_info.CellsNumDOFs[2].resize(numCell2Ds);
    mesh_dof_info.CellsBoundaryInfo[2].resize(numCell2Ds);

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            if (mesh.Cell2DNumberVertices(c) == 3)
                mesh_dof_info.CellsNumDOFs[2][c] =
                    reference_element_data.FEM_ReferenceElement_Data_3D.tetrahedron_reference_element_data.NumDofs2D;
            else if (mesh.Cell2DNumberVertices(c) == 4)
                mesh_dof_info.CellsNumDOFs[2][c] =
                    reference_element_data.FEM_ReferenceElement_Data_3D.hexahedron_reference_element_data.NumDofs2D;
            else
                throw std::runtime_error("not valid element");
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[2][c] = reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs2D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
        mesh_dof_info.CellsBoundaryInfo[2][c] = boundary_info.at(mesh.Cell2DMarker(c));
    }

    const unsigned int numCell3Ds = mesh.Cell3DTotalNumber();

    mesh_dof_info.CellsNumDOFs[3].resize(numCell3Ds);
    mesh_dof_info.CellsBoundaryInfo[3].resize(numCell3Ds);

    for (unsigned int c = 0; c < numCell3Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            if (mesh.Cell3DNumberVertices(c) == 4 && mesh.Cell3DNumberEdges(c) == 6 && mesh.Cell3DNumberFaces(c) == 4)
                mesh_dof_info.CellsNumDOFs[3][c] =
                    reference_element_data.FEM_ReferenceElement_Data_3D.tetrahedron_reference_element_data.NumDofs3D;
            else if (mesh.Cell3DNumberVertices(c) == 8 && mesh.Cell3DNumberEdges(c) == 12 && mesh.Cell3DNumberFaces(c) == 6)
                mesh_dof_info.CellsNumDOFs[3][c] =
                    reference_element_data.FEM_ReferenceElement_Data_3D.hexahedron_reference_element_data.NumDofs3D;
            else
                throw std::runtime_error("not valid element");
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[3][c] = reference_element_data.VEM_ReferenceElement_Data_3D.NumDofs3D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
        mesh_dof_info.CellsBoundaryInfo[3][c] = boundary_info.at(mesh.Cell3DMarker(c));
    }

    return mesh_dof_info;
}
//***************************************************************************
Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    Performance_Data performance;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.FEM_LocalSpace_Data.InternalQuadrature.Weights.size();
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.VEM_Performance_Data.Analysis =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_3D(),
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
} // namespace LocalSpace_PCC_3D
} // namespace PDETools
} // namespace Polydim
