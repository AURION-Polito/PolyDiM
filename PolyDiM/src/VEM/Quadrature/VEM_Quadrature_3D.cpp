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

#include "VEM_Quadrature_3D.hpp"

#include "MapTetrahedron.hpp"

#include "Quadrature_Gauss2D_Triangle.hpp"
#include "Quadrature_Gauss3D_Tetrahedron_PositiveWeights.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Quadrature
{
//****************************************************************************
VEM_QuadratureData_3D VEM_Quadrature_3D::Compute_PCC_3D(const unsigned int order) const
{
    VEM_QuadratureData_3D data;
    data.ReferenceTetrahedronQuadrature =
        Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);
    return data;
}
//****************************************************************************
VEM_QuadratureData_3D VEM_Quadrature_3D::Compute_DF_PCC_3D(const unsigned int order) const
{
    VEM_QuadratureData_3D data;
    data.ReferenceTetrahedronQuadrature =
        Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);
    return data;
}
//****************************************************************************
VEM_QuadratureData_3D VEM_Quadrature_3D::Compute_MCC_3D(const unsigned int order) const
{
    VEM_QuadratureData_3D data;
    data.ReferenceTetrahedronQuadrature =
        Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * (order + 1));
    data.QuadratureData_2D.ReferenceTriangleQuadrature =
        Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order);
    return data;
}
//****************************************************************************
Gedim::Quadrature::QuadratureData VEM_Quadrature_3D::PolyhedronInternalQuadrature(const VEM_QuadratureData_3D &data,
                                                                                  const Gedim::GeometryUtilities &geometryUtility,
                                                                                  const std::vector<Eigen::MatrixXd> &polyhedronTetrahedronVertices) const
{
    Gedim::Quadrature::QuadratureData result;

    const unsigned int numPolyhedronTetrahedrons = polyhedronTetrahedronVertices.size();

    const unsigned int numTetrahedronQuadraturePoints = data.ReferenceTetrahedronQuadrature.Points.cols();
    const unsigned int numQuadraturePoints = numPolyhedronTetrahedrons * numTetrahedronQuadraturePoints;

    result.Points.setZero(3, numQuadraturePoints);
    result.Weights.setZero(numQuadraturePoints);

    Gedim::MapTetrahedron mapTetrahedron(geometryUtility);

    for (unsigned int t = 0; t < numPolyhedronTetrahedrons; t++)
    {
        const Eigen::MatrixXd &tetrahedronVertices = polyhedronTetrahedronVertices[t];

        Gedim::MapTetrahedron::MapTetrahedronData mapTetrahedronData = mapTetrahedron.Compute(tetrahedronVertices);
        result.Points.block(0, numTetrahedronQuadraturePoints * t, 3, numTetrahedronQuadraturePoints) =
            mapTetrahedron.F(mapTetrahedronData, data.ReferenceTetrahedronQuadrature.Points);
        result.Weights.segment(numTetrahedronQuadraturePoints * t, numTetrahedronQuadraturePoints) =
            data.ReferenceTetrahedronQuadrature.Weights.array() *
            mapTetrahedron.DetJ(mapTetrahedronData, data.ReferenceTetrahedronQuadrature.Points).array().abs();
    }

    return result;
}
//****************************************************************************
VEM_Quadrature_3D::Faces_QuadratureData_PCC VEM_Quadrature_3D::PolyhedronFacesQuadrature(
    const Gedim::GeometryUtilities &geometryUtility,
    const std::vector<Eigen::MatrixXi> &polyhedronFaces,
    const std::vector<Eigen::Matrix3d> &facesRotationMatrix,
    const std::vector<Eigen::Vector3d> &facesTranslation,
    const std::vector<Eigen::Vector3d> &facesNormals,
    const std::vector<bool> &faceNormalDirections,
    const std::vector<Eigen::MatrixXd> &facesQuadraturePoints,
    const std::vector<Eigen::VectorXd> &facesQuadratureWeights) const
{
    Faces_QuadratureData_PCC result;

    const unsigned int numFaces = polyhedronFaces.size();

    unsigned int numQuadraturePoints = 0;
    for (unsigned int f = 0; f < numFaces; f++)
        numQuadraturePoints += facesQuadraturePoints[f].cols();

    result.Quadrature.Points.setZero(3, numQuadraturePoints);
    result.Quadrature.Weights.setZero(numQuadraturePoints);
    result.WeightsTimesNormal.assign(3, VectorXd::Zero(numQuadraturePoints));

    unsigned int quadraturePointOffset = 0;
    for (unsigned int f = 0; f < numFaces; f++)
    {
        const unsigned int numFaceQuadraturePoints = facesQuadraturePoints[f].cols();

        result.Quadrature.Points.block(0, quadraturePointOffset, 3, numFaceQuadraturePoints) =
            geometryUtility.RotatePointsFrom2DTo3D(facesQuadraturePoints[f], facesRotationMatrix[f], facesTranslation[f]);
        result.Quadrature.Weights.segment(quadraturePointOffset, numFaceQuadraturePoints) = facesQuadratureWeights[f];

        const double faceNormalDirection = faceNormalDirections[f] ? 1.0 : -1.0;

        for (unsigned int d = 0; d < 3; d++)
            result.WeightsTimesNormal[d].segment(quadraturePointOffset, numFaceQuadraturePoints) =
                facesQuadratureWeights[f] * faceNormalDirection * facesNormals[f](d);

        quadraturePointOffset += numFaceQuadraturePoints;
    }

    return result;
}
//****************************************************************************
VEM_Quadrature_3D::Faces_QuadratureData_MCC VEM_Quadrature_3D::PolyhedronFacesQuadrature(
    const VEM_QuadratureData_3D &data,
    const Gedim::GeometryUtilities &geometryUtility,
    const std::vector<std::vector<Eigen::Matrix3d>> &facesTriangulations2D,
    const std::vector<Eigen::Matrix3d> &facesRotationMatrix,
    const std::vector<Eigen::Vector3d> &facesTranslation) const
{
    Faces_QuadratureData_MCC result;
    VEM_Quadrature_2D quadrature2D;

    const unsigned int numFaces = facesTriangulations2D.size();
    result.FacesQuadrature.resize(numFaces);

    unsigned int numQuadraturePoints = 0;
    for (unsigned int f = 0; f < numFaces; f++)
    {
        result.FacesQuadrature[f] = quadrature2D.PolygonInternalQuadrature(data.QuadratureData_2D.ReferenceTriangleQuadrature,
                                                                           facesTriangulations2D[f]);

        numQuadraturePoints += result.FacesQuadrature[f].Points.cols();
    }

    result.Quadrature.Points.setZero(3, numQuadraturePoints);
    result.Quadrature.Weights.setZero(numQuadraturePoints);

    unsigned int quadraturePointOffset = 0;
    for (unsigned int f = 0; f < numFaces; f++)
    {
        const unsigned int numFaceQuadraturePoints = result.FacesQuadrature[f].Points.cols();

        result.Quadrature.Points.block(0, quadraturePointOffset, 3, numFaceQuadraturePoints) =
            geometryUtility.RotatePointsFrom2DTo3D(result.FacesQuadrature[f].Points, facesRotationMatrix[f], facesTranslation[f]);
        result.Quadrature.Weights.segment(quadraturePointOffset, numFaceQuadraturePoints) = result.FacesQuadrature[f].Weights;

        quadraturePointOffset += numFaceQuadraturePoints;
    }

    return result;
}
//****************************************************************************
Eigen::MatrixXd VEM_Quadrature_3D::PolyhedronInternalEdgesQuadraturePoints(const Eigen::MatrixXd &referenceSegmentInternalPoints,
                                                                           const Eigen::MatrixXd &polyhedronVertices,
                                                                           const Eigen::MatrixXi &polyhedronEdges,
                                                                           const std::vector<bool> &edgeDirections,
                                                                           const MatrixXd &edgeTangents) const
{
    const unsigned int numEdges = polyhedronEdges.cols();

    const unsigned int numEdgeInternalQuadraturePoints = referenceSegmentInternalPoints.cols();
    Eigen::MatrixXd points = MatrixXd::Zero(3, numEdgeInternalQuadraturePoints * numEdges);

    for (unsigned int i = 0; i < numEdges; i++)
    {
        const Vector3d &edgeStart = edgeDirections[i] ? polyhedronVertices.col(polyhedronEdges(0, i))
                                                      : polyhedronVertices.col(polyhedronEdges(1, i));

        const Vector3d &edgeTangent = edgeTangents.col(i);
        const double direction = edgeDirections[i] ? 1.0 : -1.0;

        for (unsigned int q = 0; q < numEdgeInternalQuadraturePoints; q++)
        {
            points.col(q + i * numEdgeInternalQuadraturePoints) =
                edgeStart + direction * referenceSegmentInternalPoints(0, q) * edgeTangent;
        }
    }

    return points;
}
//****************************************************************************
} // namespace Quadrature
} // namespace VEM
} // namespace Polydim
