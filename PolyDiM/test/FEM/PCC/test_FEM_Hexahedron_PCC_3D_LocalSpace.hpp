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

#ifndef __TEST_FEM_Hexahedron_PCC_3D_LocalSpace_H
#define __TEST_FEM_Hexahedron_PCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Hexahedron_PCC_3D_LocalSpace.hpp"
#include "FEM_Hexahedron_PCC_3D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_FEM_PCC_3D_Hexahedron_Geometry final
{
    struct Face_2D_Geometry final
    {
        Eigen::MatrixXd Vertices;
        std::vector<bool> EdgesDirection;
        Eigen::MatrixXd EdgesTangent;
        Eigen::VectorXd EdgesLength;
        std::vector<Eigen::Matrix3d> TriangulationVertices;
    };

    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    std::vector<bool> EdgesDirection;
    std::vector<double> FacesArea;
    std::vector<bool> FacesDirection;
    std::vector<Face_2D_Geometry> Faces_2D_Geometry;
    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
    std::vector<Eigen::Vector3d> FacesNormal;
    std::vector<bool> FacesNormalDirection;
    std::vector<Eigen::MatrixXd> TetrahedronVertices;
    double Volume;
};

Test_FEM_PCC_3D_Hexahedron_Geometry Test_FEM_PCC_3D_Hexa_Geometry(const Gedim::GeometryUtilities &geometry_utilities,
                                                                  const unsigned int hexahedron_type)
{
    Test_FEM_PCC_3D_Hexahedron_Geometry result;

    switch (hexahedron_type)
    {
    case 0: // reference tetrahedron
    {
        result.Vertices.resize(3, 8);
        result.Vertices << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        result.Edges.resize(2, 12);
        result.Edges.col(0) << 0, 1;
        result.Edges.col(1) << 1, 2;
        result.Edges.col(2) << 3, 2;
        result.Edges.col(3) << 0, 3;
        result.Edges.col(4) << 4, 5;
        result.Edges.col(5) << 5, 6;
        result.Edges.col(6) << 7, 6;
        result.Edges.col(7) << 4, 7;
        result.Edges.col(8) << 0, 4;
        result.Edges.col(9) << 1, 5;
        result.Edges.col(10) << 3, 7;
        result.Edges.col(11) << 2, 6;

        result.Faces.resize(6, Eigen::MatrixXi(2, 4));
        result.Faces[0].row(0) << 0, 1, 2, 3;
        result.Faces[0].row(1) << 0, 1, 2, 3;
        result.Faces[1].row(0) << 4, 5, 6, 7;
        result.Faces[1].row(1) << 4, 5, 6, 7;
        result.Faces[2].row(0) << 0, 1, 5, 4;
        result.Faces[2].row(1) << 0, 9, 4, 8;
        result.Faces[3].row(0) << 3, 2, 6, 7;
        result.Faces[3].row(1) << 2, 11, 6, 10;
        result.Faces[4].row(0) << 0, 3, 7, 4;
        result.Faces[4].row(1) << 3, 10, 7, 8;
        result.Faces[5].row(0) << 1, 2, 6, 5;
        result.Faces[5].row(1) << 1, 11, 5, 9;
    }
    break;
    case 1: {
        const Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
        const Eigen::Vector3d v2 = Eigen::Vector3d(4.0, 0.0, 0.0);
        const Eigen::Vector3d v3 = Eigen::Vector3d(0.0, 0.0, 0.5);
        const Eigen::Vector3d v4 = Eigen::Vector3d(0.0, 3.0, 0.0);
        Gedim::GeometryUtilities::Polyhedron parallelepiped = geometry_utilities.CreateParallelepipedWithOrigin(v1, v2, v3, v4);

        result.Vertices = parallelepiped.Vertices;
        result.Edges = parallelepiped.Edges;
        result.Faces = parallelepiped.Faces;
    }
    break;
    default:
        throw std::runtime_error("unknown tetrahedron type");
    }

    result.EdgesDirection.resize(12, true);
    result.FacesDirection.resize(6, true);

    Gedim::MeshUtilities mesh_utilities;

    const std::vector<unsigned int> vertexMarkers(8, 1);
    const std::vector<unsigned int> edgeMarkers(12, 1);
    const std::vector<unsigned int> faceMarkers(6, 1);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.Mesh3DFromPolyhedron(result.Vertices, result.Edges, result.Faces, vertexMarkers, edgeMarkers, faceMarkers, mesh);
    mesh_utilities.ComputeCell2DCell3DNeighbours(mesh);
    const auto geometric_data = mesh_utilities.FillMesh3DGeometricData(geometry_utilities, mesh);

    result.Faces_2D_Geometry.resize(6);
    for (unsigned int f = 0; f < 6; f++)
    {
        result.Faces_2D_Geometry[f].Vertices = geometric_data.Cell3DsFaces2DVertices[0][f];
        result.Faces_2D_Geometry[f].EdgesDirection = geometric_data.Cell3DsFacesEdgeDirections[0][f];
        result.Faces_2D_Geometry[f].EdgesTangent = geometric_data.Cell3DsFacesEdge2DTangents[0][f];
        result.Faces_2D_Geometry[f].EdgesLength = geometric_data.Cell3DsFacesEdgeLengths[0][f];
        result.Faces_2D_Geometry[f].TriangulationVertices = geometric_data.Cell3DsFaces2DTriangulations[0][f];
    }

    result.Vertices = geometric_data.Cell3DsVertices[0];
    result.Edges = geometric_data.Cell3DsEdges[0];
    result.Faces = geometric_data.Cell3DsFaces[0];
    result.EdgesDirection = geometric_data.Cell3DsEdgeDirections[0];
    result.FacesArea = geometric_data.Cell3DsFacesAreas[0];
    result.FacesDirection = geometric_data.Cell3DsFacesNormalGlobalDirection[0];
    result.FacesRotationMatrix = geometric_data.Cell3DsFacesRotationMatrices[0];
    result.FacesTranslation = geometric_data.Cell3DsFacesTranslations[0];
    result.FacesNormal = geometric_data.Cell3DsFacesNormals[0];
    result.FacesNormalDirection = geometric_data.Cell3DsFacesNormalDirections[0];
    result.TetrahedronVertices = geometric_data.Cell3DsTetrahedronPoints[0];
    result.Volume = geometric_data.Cell3DsVolumes[0];

    return result;
}

TEST(Test_FEM_Hexahedron_PCC_3D, Test_FEM_Hexahedron_PCC_3D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_Hexahedron_PCC_3D_ReferenceElement reference_element;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Test_FEM_PCC_3D_Hexahedron_Geometry geom_data = Test_FEM_PCC_3D_Hexa_Geometry(geometry_utilities, 0);

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss3D_Hexahedron::FillPointsAndWeights(10);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    for (unsigned int k = 1; k < 4; k++)
    {
        const auto reference_element_data = reference_element.Create(k);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3, dofs.cols() + referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points, reference_element_data);

        ASSERT_TRUE((dofs.leftCols(8) - reference_element.Vertices).norm() < 1.0e-13);
        ASSERT_TRUE((basisValues.topRows(dofs.cols()) - Eigen::MatrixXd::Identity(dofs.cols(), dofs.cols())).norm() < 1.0e-13);

        const std::vector<Eigen::MatrixXd> gradBasisValues =
            reference_element.EvaluateBasisFunctionDerivatives(points, reference_element_data);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumGradZValues = gradBasisValues[2].rowwise().sum();
        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE((basisValues.topRows(reference_element_data.NumBasisFunctions) -
                         Eigen::MatrixXd::Identity(reference_element_data.NumBasisFunctions, reference_element_data.NumBasisFunctions))
                            .norm() < 1.0e-13);
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradZValues[q]) < 1.0e-13);
        }
    }
}

TEST(Test_FEM_Hexahedron_PCC_3D, Test_FEM_Hexahedron_PCC_3D)
{
    const Polydim::FEM::PCC::FEM_Hexahedron_PCC_3D_ReferenceElement reference_element;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_FEM_PCC_3D_Hexahedron_Geometry hexa_data = Test_FEM_PCC_3D_Hexa_Geometry(geometry_utilities, 1);

    Polydim::FEM::PCC::FEM_PCC_3D_Polyhedron_Geometry hexa_geometry = {geometry_utilities_config.Tolerance1D,
                                                                       geometry_utilities_config.Tolerance2D,
                                                                       geometry_utilities_config.Tolerance3D,
                                                                       hexa_data.Vertices,
                                                                       hexa_data.Edges,
                                                                       hexa_data.Faces,
                                                                       {},
                                                                       hexa_data.EdgesDirection,
                                                                       hexa_data.FacesDirection,
                                                                       hexa_data.FacesRotationMatrix,
                                                                       hexa_data.FacesTranslation};
    hexa_geometry.Faces_2D_Geometry.resize(6);
    for (unsigned int f = 0; f < 6; ++f)
    {
        const auto &face_geometry = hexa_data.Faces_2D_Geometry[f];
        hexa_geometry.Faces_2D_Geometry[f] = {geometry_utilities_config.Tolerance1D,
                                              geometry_utilities_config.Tolerance2D,
                                              face_geometry.Vertices,
                                              face_geometry.EdgesDirection,
                                              face_geometry.EdgesTangent,
                                              face_geometry.EdgesLength};
    }

    const auto polyedron_faces_normal = hexa_data.FacesNormal;
    const auto polyhedron_faces_normal_direction = hexa_data.FacesNormalDirection;

    for (unsigned int k = 1; k < 4; k++)
    {
        const Polydim::FEM::PCC::FEM_Hexahedron_PCC_3D_ReferenceElement reference_element;
        const auto reference_element_data = reference_element.Create(k);

        Polydim::FEM::PCC::FEM_Hexahedron_PCC_3D_LocalSpace local_space;
        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, hexa_geometry);

        const auto &internal_quadrature = local_space_data.InternalQuadrature;

        const std::vector<Eigen::MatrixXd> gradBasisValues_ref =
            reference_element.EvaluateBasisFunctionDerivatives(local_space_data.InternalQuadrature.Points, reference_element_data);

        const Eigen::MatrixXd basis_function_values =
            local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_data);
        const std::vector<Eigen::MatrixXd> gradBasisValues =
            local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const Eigen::VectorXd sumBasisValues = basis_function_values.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumGradZValues = gradBasisValues[2].rowwise().sum();
        for (unsigned int q = 0; q < internal_quadrature.Points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradZValues[q]) < 1.0e-13);
        }

        const auto &derivative_values = local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data,
                                                                                          local_space_data,
                                                                                          internal_quadrature.Points);

        Eigen::VectorXd internal_integral = Eigen::VectorXd::Zero(reference_element_data.NumBasisFunctions);
        for (unsigned int dim = 0; dim < reference_element_data.Dimension; ++dim)
            internal_integral += derivative_values[dim].transpose() * internal_quadrature.Weights;

        if (k == 1)
        {
            Eigen::VectorXd exact_internal_integral = Eigen::VectorXd::Zero(local_space_data.NumberOfBasisFunctions);
            exact_internal_integral[0] = -3.875;
            exact_internal_integral[1] = -3.125;
            exact_internal_integral[2] = -2.125;
            exact_internal_integral[3] = -2.875;
            exact_internal_integral[4] = 2.125;
            exact_internal_integral[5] = 2.875;
            exact_internal_integral[6] = 3.875;
            exact_internal_integral[7] = 3.125;

            ASSERT_TRUE((internal_integral - exact_internal_integral).norm() <
                        1.0e-14 * std::max(1.0, exact_internal_integral.norm()));
        }

        Eigen::VectorXd boundary_integral = Eigen::VectorXd::Zero(local_space_data.NumberOfBasisFunctions);
        for (unsigned int b = 0; b < polyedron_faces_normal.size(); ++b)
        {
            const Eigen::Vector3d boundary_normal =
                polyhedron_faces_normal_direction[b] ? +1.0 * polyedron_faces_normal[b] : -1.0 * polyedron_faces_normal[b];
            const auto &boundary_quadrature = local_space_data.BoundaryQuadrature[b];
            const auto boundary_values =
                local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_data, boundary_quadrature.Points);
            boundary_integral += boundary_values.transpose() * boundary_quadrature.Weights * boundary_normal.sum();
        }

        Eigen::VectorXd error = internal_integral - boundary_integral;
        ASSERT_TRUE((internal_integral - boundary_integral).norm() < 1.0e-14 * std::max(1.0, boundary_integral.norm()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
