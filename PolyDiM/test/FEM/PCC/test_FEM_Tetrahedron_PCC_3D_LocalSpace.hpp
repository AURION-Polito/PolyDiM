#ifndef __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H
#define __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"
#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"
#include "VEM_Monomials_3D.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_FEM_PCC_3D_Tetrahedron_Geometry final
{
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    std::vector<bool> EdgesDirection;
    std::vector<double> FacesArea;
    std::vector<bool> FacesDirection;
    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
};

Test_FEM_PCC_3D_Tetrahedron_Geometry Test_FEM_PCC_3D_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_FEM_PCC_3D_Tetrahedron_Geometry result;

    const Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d v2 = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d v3 = Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d v4 = Eigen::Vector3d(0.0, 1.0, 0.0);
    Gedim::GeometryUtilities::Polyhedron tetrahedron = geometry_utilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

    result.Vertices = tetrahedron.Vertices;
    result.Edges = tetrahedron.Edges;
    result.Faces = tetrahedron.Faces;

    result.EdgesDirection.resize(6, true);
    result.FacesDirection.resize(4, true);

    const std::vector<Eigen::MatrixXd> facesVertices = geometry_utilities.PolyhedronFaceVertices(result.Vertices, result.Faces);
    result.FacesTranslation = geometry_utilities.PolyhedronFaceTranslations(facesVertices);
    const auto FacesNormals = geometry_utilities.PolyhedronFaceNormals(facesVertices);
    result.FacesRotationMatrix =
        geometry_utilities.PolyhedronFaceRotationMatrices(facesVertices, FacesNormals, result.FacesTranslation);

    const unsigned int numFaces = result.Faces.size();
    result.FacesArea.resize(numFaces);

    for (unsigned int f = 0; f < numFaces; f++)
    {
        const auto faceVertices = geometry_utilities.RotatePointsFrom3DTo2D(facesVertices[f],
                                                                            result.FacesRotationMatrix[f].transpose(),
                                                                            result.FacesTranslation[f]);
        result.FacesArea[f] = geometry_utilities.PolygonArea(faceVertices);
    }

    return result;
}

TEST(Test_FEM_Tetrahedron_PCC_3D, Test_FEM_Tetrahedron_PCC_3D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

    for (unsigned int k = 1; k < 4; k++)
    {

        const auto referenceQuadrature =
            Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(10);
        const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

        const auto reference_element_data = reference_element.Create(k);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3, dofs.cols() + referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points, reference_element_data);
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

        VEM::Monomials::VEM_Monomials_3D monomials;
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

TEST(Test_FEM_Tetrahedron_PCC_3D, Test_FEM_Tetrahedron_PCC_3D)
{
    const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_FEM_PCC_3D_Tetrahedron_Geometry tetrahedron_data = Test_FEM_PCC_3D_Geometry(geometry_utilities);

    Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_Geometry tetra_geometry = {geometry_utilities_config.Tolerance1D,
                                                                         geometry_utilities_config.Tolerance2D,
                                                                         geometry_utilities_config.Tolerance3D,
                                                                         tetrahedron_data.Vertices,
                                                                         tetrahedron_data.Edges,
                                                                         tetrahedron_data.Faces,
                                                                         tetrahedron_data.EdgesDirection,
                                                                         tetrahedron_data.FacesArea,
                                                                         tetrahedron_data.FacesDirection,
                                                                         tetrahedron_data.FacesRotationMatrix,
                                                                         tetrahedron_data.FacesTranslation};

    for (unsigned int k = 1; k < 4; k++)
    {
        const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;
        const auto reference_element_data = reference_element.Create(k);

        Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_LocalSpace local_space;
        const auto local_space_Data = local_space.CreateLocalSpace(reference_element_data, tetra_geometry);

        const auto quadrature = local_space_Data.InternalQuadrature;

        const std::vector<Eigen::MatrixXd> gradBasisValues_ref =
            reference_element.EvaluateBasisFunctionDerivatives(local_space_Data.InternalQuadrature.Points, reference_element_data);

        const Eigen::MatrixXd basis_function_values =
            local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_Data);
        const std::vector<Eigen::MatrixXd> gradBasisValues =
            local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_Data);

        const Eigen::VectorXd sumBasisValues = basis_function_values.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumGradZValues = gradBasisValues[2].rowwise().sum();
        for (unsigned int q = 0; q < quadrature.Points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-13);
            ASSERT_TRUE(abs(sumGradZValues[q]) < 1.0e-13);
        }

        VEM::Monomials::VEM_Monomials_3D monomials;
        const auto monomials_data = monomials.Compute(k);
        const Eigen::MatrixXd vander_matrix = monomials.Vander(monomials_data, quadrature.Points, Eigen::Vector3d::Zero(), 1.0);

        const std::vector<Eigen::MatrixXd> grad_vander_matrix = monomials.VanderDerivatives(monomials_data, vander_matrix, 1.0);

        const Eigen::MatrixXd Hmatrix = vander_matrix.transpose() * quadrature.Weights.asDiagonal() * vander_matrix;
        const Eigen::MatrixXd rhs = vander_matrix.transpose() * quadrature.Weights.asDiagonal() * basis_function_values;
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
