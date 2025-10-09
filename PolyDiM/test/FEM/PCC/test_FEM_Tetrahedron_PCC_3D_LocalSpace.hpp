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

#ifndef __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H
#define __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"
#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"
#include "Monomials_3D.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_FEM_PCC_3D_Tetrahedron_Geometry final
{
    struct Face_2D_Geometry final
    {
        Eigen::MatrixXd Vertices;
        std::vector<bool> EdgesDirection;
        Eigen::MatrixXd EdgesTangent;
        Eigen::VectorXd EdgesLength;
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
    double Volume;
};

Test_FEM_PCC_3D_Tetrahedron_Geometry Test_FEM_PCC_3D_Geometry(const Gedim::GeometryUtilities &geometry_utilities,
                                                              const unsigned int tetrahedron_type)
{
    Test_FEM_PCC_3D_Tetrahedron_Geometry result;

    switch (tetrahedron_type)
    {
    case 0: // reference tetrahedron
    {
        result.Vertices.resize(3, 4);
        result.Vertices.col(0) << 0.0, 0.0, 0.0;
        result.Vertices.col(1) << 1.0, 0.0, 0.0;
        result.Vertices.col(2) << 0.0, 1.0, 0.0;
        result.Vertices.col(3) << 0.0, 0.0, 1.0;

        result.Edges.resize(2, 6);
        result.Edges.col(0) << 0, 1;
        result.Edges.col(1) << 1, 2;
        result.Edges.col(2) << 2, 0;
        result.Edges.col(3) << 0, 3;
        result.Edges.col(4) << 1, 3;
        result.Edges.col(5) << 2, 3;

        result.Faces.resize(4, Eigen::MatrixXi(2, 3));
        result.Faces[0].row(0) << 0, 1, 2;
        result.Faces[0].row(1) << 0, 1, 2;
        result.Faces[1].row(0) << 0, 1, 3;
        result.Faces[1].row(1) << 0, 4, 3;
        result.Faces[2].row(0) << 0, 2, 3;
        result.Faces[2].row(1) << 2, 5, 3;
        result.Faces[3].row(0) << 1, 2, 3;
        result.Faces[3].row(1) << 1, 5, 4;
    }
    break;
    case 1: {
        const Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
        const Eigen::Vector3d v2 = Eigen::Vector3d(1.0, 0.0, 0.0);
        const Eigen::Vector3d v3 = Eigen::Vector3d(0.0, 0.0, 1.0);
        const Eigen::Vector3d v4 = Eigen::Vector3d(0.0, 1.0, 0.0);
        Gedim::GeometryUtilities::Polyhedron tetrahedron = geometry_utilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

        result.Vertices = tetrahedron.Vertices;
        result.Edges = tetrahedron.Edges;
        result.Faces = tetrahedron.Faces;
    }
    break;
    default:
        throw std::runtime_error("unknown tetrahedron type");
    }

    result.EdgesDirection.resize(6, true);
    result.FacesDirection.resize(4, true);

    const std::vector<Eigen::MatrixXd> facesVertices = geometry_utilities.PolyhedronFaceVertices(result.Vertices, result.Faces);
    result.FacesTranslation = geometry_utilities.PolyhedronFaceTranslations(facesVertices);
    result.FacesNormal = geometry_utilities.PolyhedronFaceNormals(facesVertices);

    const Eigen::Vector3d barycenter = geometry_utilities.PolyhedronBarycenter(result.Vertices);
    result.FacesNormalDirection = geometry_utilities.PolyhedronFaceNormalDirections(facesVertices, barycenter, result.FacesNormal);
    result.FacesRotationMatrix =
        geometry_utilities.PolyhedronFaceRotationMatrices(facesVertices, result.FacesNormal, result.FacesTranslation);

    const auto faces_rotated_vertices =
        geometry_utilities.PolyhedronFaceRotatedVertices(facesVertices, result.FacesTranslation, result.FacesRotationMatrix);
    const unsigned int numFaces = result.Faces.size();
    result.FacesArea.resize(numFaces);
    result.Faces_2D_Geometry.resize(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
    {
        auto &face_2D_geometry = result.Faces_2D_Geometry[f];

        face_2D_geometry.Vertices = faces_rotated_vertices[f];
        face_2D_geometry.EdgesDirection.resize(face_2D_geometry.Vertices.cols(), true);
        face_2D_geometry.EdgesLength = geometry_utilities.PolygonEdgeLengths(face_2D_geometry.Vertices);
        face_2D_geometry.EdgesTangent = geometry_utilities.PolygonEdgeTangents(face_2D_geometry.Vertices);
        result.FacesArea[f] = geometry_utilities.PolygonArea(face_2D_geometry.Vertices);
    }

    result.Volume = geometry_utilities.PolyhedronVolumeByInternalIntegral({result.Vertices});

    return result;
}

TEST(Test_FEM_Tetrahedron_PCC_3D, Test_FEM_Tetrahedron_PCC_3D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(10);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    for (unsigned int k = 1; k < 6; k++)
    {
        const auto reference_element_data = reference_element.Create(k);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3, dofs.cols() + referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points, reference_element_data);

        ASSERT_TRUE((basisValues.topRows(reference_element_data.NumBasisFunctions) -
                     Eigen::MatrixXd::Identity(reference_element_data.NumBasisFunctions, reference_element_data.NumBasisFunctions))
                        .norm() < 1.0e-13);

        const std::vector<Eigen::MatrixXd> gradBasisValues =
            reference_element.EvaluateBasisFunctionDerivatives(referenceQuadraturePoints, reference_element_data);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumGradZValues = gradBasisValues[2].rowwise().sum();
        for (unsigned int q = 0; q < referenceQuadraturePoints.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q + dofs.cols()] - 1.0) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradZValues[q]) < 1.0e-13);
        }

        Utilities::Monomials_3D monomials;
        const auto monomials_data = monomials.Compute(k);
        const Eigen::MatrixXd vander_matrix =
            monomials.Vander(monomials_data, referenceQuadraturePoints, Eigen::Vector3d::Zero(), 1.0);

        const std::vector<Eigen::MatrixXd> grad_vander_matrix = monomials.VanderDerivatives(monomials_data, vander_matrix, 1.0);

        const Eigen::MatrixXd Hmatrix = vander_matrix.transpose() * referenceQuadrature.Weights.asDiagonal() * vander_matrix;
        const Eigen::MatrixXd rhs = vander_matrix.transpose() * referenceQuadrature.Weights.asDiagonal() *
                                    basisValues.bottomRows(referenceQuadraturePoints.cols());
        const Eigen::MatrixXd coefficients = Hmatrix.llt().solve(rhs);

        for (unsigned int d = 0; d < 3; d++)
        {
            const double norm_der = gradBasisValues[d].bottomRows(referenceQuadraturePoints.cols()).norm();
            ASSERT_TRUE(
                (gradBasisValues[d].bottomRows(referenceQuadraturePoints.cols()) - grad_vander_matrix[d] * coefficients).norm() <
                1.0e-9 * norm_der);
        }
    }
}

TEST(Test_FEM_Tetrahedron_PCC_3D, Test_FEM_Tetrahedron_PCC_3D_local_space)
{
    const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_FEM_PCC_3D_Tetrahedron_Geometry tetrahedron_data = Test_FEM_PCC_3D_Geometry(geometry_utilities, 1);

    Polydim::FEM::PCC::FEM_PCC_3D_Polyhedron_Geometry tetra_geometry = {geometry_utilities_config.Tolerance1D,
                                                                        geometry_utilities_config.Tolerance2D,
                                                                        geometry_utilities_config.Tolerance3D,
                                                                        tetrahedron_data.Vertices,
                                                                        tetrahedron_data.Edges,
                                                                        tetrahedron_data.Faces,
                                                                        tetrahedron_data.EdgesDirection,
                                                                        tetrahedron_data.FacesDirection,
                                                                        tetrahedron_data.FacesRotationMatrix,
                                                                        tetrahedron_data.FacesTranslation,
                                                                        {}};
    tetra_geometry.Faces_2D_Geometry.resize(4);
    for (unsigned int f = 0; f < 4; ++f)
    {
        const auto &face_geometry = tetrahedron_data.Faces_2D_Geometry[f];
        tetra_geometry.Faces_2D_Geometry[f] = {geometry_utilities_config.Tolerance1D,
                                               geometry_utilities_config.Tolerance2D,
                                               face_geometry.Vertices,
                                               face_geometry.EdgesDirection,
                                               face_geometry.EdgesTangent,
                                               face_geometry.EdgesLength};
    }

    const auto polyedron_faces_normal = tetrahedron_data.FacesNormal;
    const auto polyhedron_faces_normal_direction = tetrahedron_data.FacesNormalDirection;
    const auto polyhedron_volume = tetrahedron_data.Volume;

    for (unsigned int k = 1; k < 5; k++)
    {
        const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;
        const auto reference_element_data = reference_element.Create(k);

        Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_LocalSpace local_space;
        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, tetra_geometry);

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
            exact_internal_integral[0] = polyhedron_volume * -3.0;
            exact_internal_integral[1] = polyhedron_volume;
            exact_internal_integral[2] = polyhedron_volume;
            exact_internal_integral[3] = polyhedron_volume;
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

        ASSERT_TRUE((internal_integral - boundary_integral).norm() < 1.0e-14 * std::max(1.0, boundary_integral.norm()));

        Utilities::Monomials_3D monomials;
        const auto monomials_data = monomials.Compute(k);
        const Eigen::MatrixXd vander_matrix =
            monomials.Vander(monomials_data, internal_quadrature.Points, Eigen::Vector3d::Zero(), 1.0);

        const std::vector<Eigen::MatrixXd> grad_vander_matrix = monomials.VanderDerivatives(monomials_data, vander_matrix, 1.0);

        const Eigen::MatrixXd Hmatrix = vander_matrix.transpose() * internal_quadrature.Weights.asDiagonal() * vander_matrix;
        const Eigen::MatrixXd rhs = vander_matrix.transpose() * internal_quadrature.Weights.asDiagonal() * basis_function_values;
        const Eigen::MatrixXd coefficients = Hmatrix.llt().solve(rhs);

        for (unsigned int d = 0; d < 3; d++)
        {
            const double norm_der = gradBasisValues[d].norm();
            ASSERT_TRUE((gradBasisValues[d] - grad_vander_matrix[d] * coefficients).norm() < 1.0e-9 * norm_der);
        }
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
