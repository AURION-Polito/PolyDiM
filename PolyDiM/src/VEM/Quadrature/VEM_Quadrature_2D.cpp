#include "VEM_Quadrature_2D.hpp"

#include "MapTriangle.hpp"

#include "Quadrature/Quadrature_Gauss2D_Triangle.hpp"
#include "Quadrature/Quadrature_GaussLobatto1D.hpp"
#include "Quadrature/Quadrature_Gauss1D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Quadrature
{
//****************************************************************************
VEM_QuadratureData_2D VEM_Quadrature_2D::Compute_PCC_2D(const unsigned int order) const
{
    VEM_QuadratureData_2D data;

    data.ReferenceTriangleQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order);

    data.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_GaussLobatto1D::FillPointsAndWeights(2 * order - 1);

    if (order == 1)
    {
        data.ReferenceSegmentInternalPoints.resize(0, 0);
        data.ReferenceSegmentInternalWeights.resize(0);
        data.ReferenceSegmentExtremaWeights.setConstant(0.5);
    }
    else
    {
        const unsigned int edgeReferenceQuadratureNumPoints = data.ReferenceSegmentQuadrature.Points.cols();

        data.ReferenceSegmentInternalPoints = data.ReferenceSegmentQuadrature.Points.block(0,
                                                                                           1,
                                                                                           3,
                                                                                           edgeReferenceQuadratureNumPoints - 2);
        data.ReferenceSegmentInternalWeights = data.ReferenceSegmentQuadrature.Weights.segment(1,
                                                                                               edgeReferenceQuadratureNumPoints - 2);
        data.ReferenceSegmentExtremaWeights<< data.ReferenceSegmentQuadrature.Weights(0),
            data.ReferenceSegmentQuadrature.Weights.tail<1>();
    }

    return data;
}
//****************************************************************************
VEM_QuadratureData_2D VEM_Quadrature_2D::Compute_MCC_2D(const unsigned int order) const
{
    VEM_QuadratureData_2D data;

    data.ReferenceTriangleQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * (order + 1));
    data.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * order + 1);
    data.ReferenceSegmentInternalPoints = data.ReferenceSegmentQuadrature.Points;
    data.ReferenceSegmentInternalWeights = data.ReferenceSegmentQuadrature.Weights;

    return data;
}
//****************************************************************************
Gedim::Quadrature::QuadratureData VEM_Quadrature_2D::PolygonInternalQuadrature(const VEM_QuadratureData_2D& data,
                                                                               const std::vector<Eigen::Matrix3d>& polygonTriangulationVertices) const
{
    Gedim::Quadrature::QuadratureData result;

    const unsigned int numPolygonTriangles = polygonTriangulationVertices.size();

    const unsigned int numTriangleQuadraturePoints = data.ReferenceTriangleQuadrature.Points.cols();
    const unsigned int numQuadraturePoints = numPolygonTriangles * numTriangleQuadraturePoints;

    result.Points.setZero(3, numQuadraturePoints);
    result.Weights.setZero(numQuadraturePoints);

    Gedim::MapTriangle mapTriangle;

    for (unsigned int t = 0; t < numPolygonTriangles; t++)
    {
        const Eigen::Matrix3d& triangleVertices = polygonTriangulationVertices[t];

        Gedim::MapTriangle::MapTriangleData mapData = mapTriangle.Compute(triangleVertices);
        result.Points.block(0,
                            numTriangleQuadraturePoints * t,
                            3,
                            numTriangleQuadraturePoints) = mapTriangle.F(mapData,
                            data.ReferenceTriangleQuadrature.Points);
        result.Weights.segment(numTriangleQuadraturePoints * t,
                               numTriangleQuadraturePoints) = data.ReferenceTriangleQuadrature.Weights.array() *
              mapTriangle.DetJ(mapData,
                               data.ReferenceTriangleQuadrature.Points).array().abs();

    }

    return result;
}
//****************************************************************************
VEM_Quadrature_2D::Edges_QuadratureData VEM_Quadrature_2D::PolygonEdgesLobattoQuadrature(const VEM_QuadratureData_2D& data,
                                                                                         const Eigen::MatrixXd& polygonVertices,
                                                                                         const Eigen::VectorXd& edgeLengths,
                                                                                         const std::vector<bool>& edgeDirections,
                                                                                         const Eigen::MatrixXd& edgeTangents,
                                                                                         const Eigen::MatrixXd& edgeNormals) const
{
    Edges_QuadratureData result;

    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    const unsigned int numVerticesQuadraturePoints = numVertices;
    const unsigned int numEdgeQuadraturePoints = data.ReferenceSegmentInternalWeights.size() * numEdges;
    const unsigned int numQuadraturePoints = numVerticesQuadraturePoints + numEdgeQuadraturePoints;

    result.Quadrature.Points.resize(3, numQuadraturePoints);
    result.Quadrature.Weights.setZero(numQuadraturePoints);
    result.WeightsTimesNormal.assign(2, VectorXd::Zero(numQuadraturePoints));

    for(unsigned int i = 0; i < numVertices; ++i)
        result.Quadrature.Points.col(i) = polygonVertices.col(i);

    // offset used below to set edge-internal quadrature points and weights.
    unsigned int edgeInternalPointsOffset = numVertices;
    for(unsigned int i = 0; i < numEdges; ++i)
    {
        const Vector2d& refEdgeExtremaWeights = data.ReferenceSegmentExtremaWeights;
        const double absMapDeterminant = std::abs(edgeLengths[i]);
        VectorXd outNormalTimesAbsMapDeterminant = edgeNormals.col(i) * absMapDeterminant;

        result.Quadrature.Weights(i) += refEdgeExtremaWeights(0) * absMapDeterminant;
        result.Quadrature.Weights((i + 1) % numEdges) += refEdgeExtremaWeights(1) * absMapDeterminant;
        for(unsigned int d = 0; d < 2; ++d)
        {
            result.WeightsTimesNormal[d](i) +=
                outNormalTimesAbsMapDeterminant(d) * refEdgeExtremaWeights(0);
            result.WeightsTimesNormal[d]((i + 1) % numEdges) +=
                outNormalTimesAbsMapDeterminant(d) * refEdgeExtremaWeights(1);
        }

        if(data.ReferenceSegmentInternalWeights.size() > 0)
        {
            // map edge internal quadrature points
            const VectorXd& refEdgeInternalWeights = data.ReferenceSegmentInternalWeights;

            const Vector3d& edgeStart = edgeDirections[i] ? polygonVertices.col(i) : polygonVertices.col((i + 1) % numVertices);
            const Vector3d& edgeTangent = edgeTangents.col(i);
            const double direction = edgeDirections[i] ? 1.0 : -1.0;

            const unsigned int numEdgeInternalQuadraturePoints = data.ReferenceSegmentInternalPoints.cols();
            MatrixXd edgeInternalQuadraturePoints(3, numEdgeInternalQuadraturePoints);
            for (unsigned int q = 0; q < numEdgeInternalQuadraturePoints; q++)
            {
                edgeInternalQuadraturePoints.col(q) = edgeStart +
                                                      direction *
                                                          data.ReferenceSegmentInternalPoints(0, q) *
                                                          edgeTangent;
            }

            result.Quadrature.Points.block(0,
                                           edgeInternalPointsOffset,
                                           3,
                                           numEdgeInternalQuadraturePoints) =
                edgeInternalQuadraturePoints;
            result.Quadrature.Weights.segment(edgeInternalPointsOffset,
                                              numEdgeInternalQuadraturePoints) =
                refEdgeInternalWeights * absMapDeterminant;
            MatrixXd edgeInternalWeightsTimesNormal =
                refEdgeInternalWeights * outNormalTimesAbsMapDeterminant.transpose();
            for(unsigned int d = 0; d < 2; ++d)
            {
                result.WeightsTimesNormal[d].segment(edgeInternalPointsOffset,
                                                     numEdgeInternalQuadraturePoints) =
                    edgeInternalWeightsTimesNormal.col(d);
            }
            edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;
        }
    }

    return result;
}
//****************************************************************************
VEM_Quadrature_2D::Edges_QuadratureData VEM_Quadrature_2D::PolygonEdgesQuadrature(const VEM_QuadratureData_2D& data,
                                                                                  const Eigen::MatrixXd& polygonVertices,
                                                                                  const Eigen::VectorXd& edgeLengths,
                                                                                  const std::vector<bool>& edgeDirections,
                                                                                  const Eigen::MatrixXd& edgeTangents,
                                                                                  const Eigen::MatrixXd& edgeNormals) const
{
    Edges_QuadratureData result;

    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    const unsigned int numQuadraturePoints = data.ReferenceSegmentInternalWeights.size() * numEdges;

    result.Quadrature.Points.resize(3, numQuadraturePoints);
    result.Quadrature.Weights.setZero(numQuadraturePoints);
    result.WeightsTimesNormal.assign(2, VectorXd::Zero(numQuadraturePoints));

    // offset used below to set edge-internal quadrature points and weights.
    unsigned int edgeInternalPointsOffset = 0;
    for(unsigned int i = 0; i < numEdges; ++i)
    {
        const double absMapDeterminant = std::abs(edgeLengths[i]);
        VectorXd outNormalTimesAbsMapDeterminant = edgeNormals.col(i) * absMapDeterminant;

        // map edge internal quadrature points
        const VectorXd& refEdgeInternalWeights = data.ReferenceSegmentInternalWeights;

        const Vector3d& edgeStart = edgeDirections[i] ? polygonVertices.col(i) : polygonVertices.col((i + 1) % numVertices);
        const Vector3d& edgeTangent = edgeTangents.col(i);
        const double direction = edgeDirections[i] ? 1.0 : -1.0;

        const unsigned int numEdgeInternalQuadraturePoints = data.ReferenceSegmentInternalPoints.cols();
        MatrixXd edgeInternalQuadraturePoints(3, numEdgeInternalQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeInternalQuadraturePoints; q++)
        {
            edgeInternalQuadraturePoints.col(q) = edgeStart +
                                                  direction *
                                                      data.ReferenceSegmentInternalPoints(0, q) *
                                                      edgeTangent;
        }

        result.Quadrature.Points.block(0,
                                       edgeInternalPointsOffset,
                                       3,
                                       numEdgeInternalQuadraturePoints) =
            edgeInternalQuadraturePoints;
        result.Quadrature.Weights.segment(edgeInternalPointsOffset,
                                          numEdgeInternalQuadraturePoints) =
            refEdgeInternalWeights * absMapDeterminant;
        MatrixXd edgeInternalWeightsTimesNormal =
            refEdgeInternalWeights * outNormalTimesAbsMapDeterminant.transpose();
        for(unsigned int d = 0; d < 2; ++d)
        {
            result.WeightsTimesNormal[d].segment(edgeInternalPointsOffset,
                                                 numEdgeInternalQuadraturePoints) =
                edgeInternalWeightsTimesNormal.col(d);
        }
        edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;

    }

    return result;
}
//****************************************************************************
}
}
}
