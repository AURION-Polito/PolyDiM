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

#include "Inertia_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace Utilities
{
// ***************************************************************************
void Inertia_Utilities::InertiaMapping2D(const MatrixXd &vertices,
                                         const Vector3d &centroid,
                                         const double &diameter,
                                         const vector<Matrix3d> &triangulation_vertices,
                                         Inertia_Data &inertia_data) const
{

    inertia_data.FmatrixInv = Matrix3d::Identity();
    inertia_data.Fmatrix = Matrix3d::Identity();
    inertia_data.absDetFmatrix = 1.0;
    inertia_data.translation = Vector3d::Zero();

    // First Rescaling
    const double invDiameter = (1.0 / diameter);

    inertia_data.FmatrixInv = invDiameter * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = inertia_data.Fmatrix * diameter;
    inertia_data.absDetFmatrix *= diameter * diameter;
    inertia_data.translation += Vector3d::Zero();

    const MatrixXd verticesFirstRescaling = inertia_data.FmatrixInv * (vertices.colwise() - inertia_data.translation);
    const Vector3d centroidFirstRescaling = inertia_data.FmatrixInv * (centroid - inertia_data.translation);

    // Inertia Mapping
    const unsigned int numTriangles = triangulation_vertices.size();
    vector<Matrix3d> triangulationsFirstRescaling(numTriangles);
    for (unsigned int n = 0; n < numTriangles; n++)
        triangulationsFirstRescaling[n] = inertia_data.FmatrixInv * (triangulation_vertices[n].colwise() - inertia_data.translation);

    const Matrix2d HmatrixFirstRescaling = geometryUtilities.PolygonMass(centroidFirstRescaling, triangulationsFirstRescaling);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(HmatrixFirstRescaling);
    if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("Eigen fails");

    const Vector2d sqrtLambdaFirstRescaling = eigensolver.eigenvalues().array().sqrt();
    const Vector2d sqrtLambdaInvFirstRescaling = eigensolver.eigenvalues().array().rsqrt();
    const MatrixXd QmatrixFirstRescaling = eigensolver.eigenvectors();
    const MatrixXd Bmatrix =
        sqrtLambdaInvFirstRescaling.asDiagonal() * QmatrixFirstRescaling.transpose() * sqrtLambdaFirstRescaling.maxCoeff();
    const MatrixXd BmatrixInv =
        QmatrixFirstRescaling * sqrtLambdaFirstRescaling.asDiagonal() / sqrtLambdaFirstRescaling.maxCoeff();

    inertia_data.FmatrixInv.topLeftCorner(2, 2) = Bmatrix * inertia_data.FmatrixInv.topLeftCorner(2, 2);
    inertia_data.Fmatrix.topLeftCorner(2, 2) = inertia_data.Fmatrix.topLeftCorner(2, 2) * BmatrixInv;
    const double detB = BmatrixInv.determinant();
    inertia_data.signDetQ = (detB > 1.0e-12) ? 1.0 : -1.0;

#ifdef TEST_INERTIA
    if (abs(detB) <= 1.0e-12)
        throw runtime_error("Singular matrix of inertia");
#endif

    inertia_data.absDetFmatrix *= abs(detB);
    inertia_data.translation += centroid;

    const MatrixXd verticesInertiaMapping = inertia_data.FmatrixInv * (vertices.colwise() - inertia_data.translation);

    // Second rescaling

    const double diameterInertiaMapping = geometryUtilities.PolygonDiameter(verticesInertiaMapping);

    const double invDiameterInertiaMapping = (1.0 / diameterInertiaMapping);

    inertia_data.FmatrixInv = invDiameterInertiaMapping * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = diameterInertiaMapping * inertia_data.Fmatrix;
    inertia_data.absDetFmatrix *= diameterInertiaMapping * diameterInertiaMapping;
    inertia_data.translation += Vector3d::Zero();
}
// ***************************************************************************
void Inertia_Utilities::InertiaMapping3D(const MatrixXd &vertices,
                                         const Vector3d &centroid,
                                         const double &diameter,
                                         const vector<MatrixXd> &tetrahedrons_vertices,
                                         Inertia_Data &inertia_data) const
{
    inertia_data.FmatrixInv = Matrix3d::Identity();
    inertia_data.Fmatrix = Matrix3d::Identity();
    inertia_data.absDetFmatrix = 1.0;
    inertia_data.translation = Vector3d::Zero();

    // First Rescaling
    const double invPolyhedronDiameter = (1.0 / diameter);

    inertia_data.FmatrixInv = invPolyhedronDiameter * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = inertia_data.Fmatrix * diameter;
    inertia_data.absDetFmatrix *= diameter * diameter * diameter;
    inertia_data.translation += Vector3d::Zero();

    const MatrixXd polyhedronVerticesFirstRescaling = inertia_data.FmatrixInv * (vertices.colwise() - inertia_data.translation);
    const Vector3d polyhedronCentroidsFirstRescaling = inertia_data.FmatrixInv * (centroid - inertia_data.translation);

    // Inertia Mapping
    const unsigned int numTetrahedron = tetrahedrons_vertices.size();
    vector<MatrixXd> polyhedronTetrahedronFirstRescaling(numTetrahedron);
    for (unsigned int n = 0; n < numTetrahedron; n++)
    {
        polyhedronTetrahedronFirstRescaling[n] =
            inertia_data.FmatrixInv * (tetrahedrons_vertices[n].colwise() - inertia_data.translation);
    }

    const Matrix3d HmatrixFirstRescaling =
        geometryUtilities.PolyhedronMass(polyhedronCentroidsFirstRescaling, polyhedronTetrahedronFirstRescaling);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(HmatrixFirstRescaling);
    if (eigensolver.info() != Eigen::Success)
        throw runtime_error("Eigen fails");

    const Vector3d sqrtLambdaFirstRescaling = eigensolver.eigenvalues().array().sqrt();
    const Vector3d sqrtLambdaInvFirstRescaling = eigensolver.eigenvalues().array().rsqrt();
    const Matrix3d QmatrixFirstRescaling = eigensolver.eigenvectors();
    const Matrix3d Bmatrix =
        sqrtLambdaInvFirstRescaling.asDiagonal() * QmatrixFirstRescaling.transpose() * sqrtLambdaFirstRescaling.maxCoeff();
    const Matrix3d BmatrixInv =
        QmatrixFirstRescaling * sqrtLambdaFirstRescaling.asDiagonal() / sqrtLambdaFirstRescaling.maxCoeff();

    inertia_data.FmatrixInv = Bmatrix * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = inertia_data.Fmatrix * BmatrixInv;
    double detB = BmatrixInv.determinant();

#ifdef TEST_INERTIA
    if (abs(detB) <= 1.0e-12)
        throw runtime_error("B singular in bulk");
#endif

    inertia_data.absDetFmatrix *= abs(detB);
    inertia_data.translation += centroid;

    const MatrixXd polyhedronVerticesInertiaMapping = inertia_data.FmatrixInv * (vertices.colwise() - inertia_data.translation);

    // Second rescaling

    const double polyhedronDiameterInertiaMapping = geometryUtilities.PolyhedronDiameter(polyhedronVerticesInertiaMapping);
    const double invPolyhedronDiameterInertiaMapping = (1.0 / polyhedronDiameterInertiaMapping);

    inertia_data.FmatrixInv = invPolyhedronDiameterInertiaMapping * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = polyhedronDiameterInertiaMapping * inertia_data.Fmatrix;
    inertia_data.absDetFmatrix *= polyhedronDiameterInertiaMapping * polyhedronDiameterInertiaMapping * polyhedronDiameterInertiaMapping;
    inertia_data.translation += Vector3d::Zero();
}
//****************************************************************************
} // namespace Utilities
} // namespace Polydim
