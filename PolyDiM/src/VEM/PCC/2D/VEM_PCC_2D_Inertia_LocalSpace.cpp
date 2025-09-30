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

#include "VEM_PCC_2D_Inertia_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Inertia_LocalSpace::CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Utilities::Inertia_Utilities inertia_utilities(geometryUtilities);

    inertia_utilities.InertiaMapping2D(polygon.Vertices,
                                       polygon.Centroid,
                                       polygon.Diameter,
                                       polygon.TriangulationVertices,
                                       localSpace.inertia_data);

    ComputeGeometryProperties(geometryUtilities, localSpace.inertia_data, polygon, localSpace.inertia_polygon);

    Quadrature::VEM_Quadrature_2D quadrature;
    Gedim::Quadrature::QuadratureData MappedInternalQuadrature =
        quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                             localSpace.inertia_polygon.TriangulationVertices);

#ifdef TEST_INERTIA
    if (abs(MappedInternalQuadrature.Weights.sum() - localSpace.inertia_data.Measure) >= 1.0e-12)
        throw runtime_error("Weights inertia are wrong - 1");
#endif

    localSpace.InternalQuadrature.Points =
        (localSpace.inertia_data.Fmatrix * MappedInternalQuadrature.Points).colwise() + localSpace.inertia_data.translation;
    localSpace.InternalQuadrature.Weights = localSpace.inertia_data.absDetFmatrix * MappedInternalQuadrature.Weights;

#ifdef TEST_INERTIA
    if (abs(localSpace.InternalQuadrature.Weights.sum() - polygon.Measure) >= 1.0e-12)
        throw runtime_error("Weights inertia are wrong - 2");
#endif

    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData MappedBoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 localSpace.inertia_polygon.Vertices,
                                                 localSpace.inertia_polygon.EdgesLength,
                                                 localSpace.inertia_polygon.EdgesDirection,
                                                 localSpace.inertia_polygon.EdgesTangent,
                                                 localSpace.inertia_polygon.EdgesNormal);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 polygon.Vertices,
                                                 polygon.EdgesLength,
                                                 polygon.EdgesDirection,
                                                 polygon.EdgesTangent,
                                                 polygon.EdgesNormal);

    InitializeProjectorsComputation(reference_element_data,
                                    localSpace.inertia_polygon.Vertices,
                                    localSpace.inertia_polygon.Centroid,
                                    localSpace.inertia_polygon.Diameter,
                                    MappedInternalQuadrature.Points,
                                    MappedInternalQuadrature.Weights,
                                    MappedBoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    ComputePiNabla(reference_element_data,
                   localSpace.inertia_polygon.Measure,
                   localSpace.inertia_polygon.Diameter,
                   MappedInternalQuadrature.Weights,
                   MappedBoundaryQuadrature.Quadrature.Weights,
                   MappedBoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(localSpace.inertia_polygon.Measure, localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data,
                                     localSpace.inertia_polygon.Measure,
                                     localSpace.inertia_polygon.Diameter,
                                     MappedBoundaryQuadrature.WeightsTimesNormal,
                                     localSpace);

    ComputePolynomialsDofs(localSpace.inertia_polygon.Measure, localSpace);

    localSpace.Diameter = localSpace.inertia_polygon.Diameter;
    localSpace.Centroid = localSpace.inertia_polygon.Centroid;
    localSpace.Measure = localSpace.inertia_polygon.Measure;

    localSpace.constantStiff = 1.0;
    localSpace.constantMass = polygon.Measure;

    return localSpace;
}
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Inertia_LocalSpace::Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                             const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Utilities::Inertia_Utilities inertia_utilities(geometryUtilities);

    inertia_utilities.InertiaMapping2D(polygon.Vertices,
                                       polygon.Centroid,
                                       polygon.Diameter,
                                       polygon.TriangulationVertices,
                                       localSpace.inertia_data);

    ComputeGeometryProperties(geometryUtilities, localSpace.inertia_data, polygon, localSpace.inertia_polygon);

    Quadrature::VEM_Quadrature_2D quadrature;
    Gedim::Quadrature::QuadratureData MappedInternalQuadrature =
        quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                             localSpace.inertia_polygon.TriangulationVertices);

#ifdef TEST
    if (abs(MappedInternalQuadrature.Weights.sum() - localSpace.inertia_data.Measure) >= 1.0e-12)
        throw runtime_error("Weights inertia are wrong - 1");
#endif

    localSpace.InternalQuadrature.Points =
        (localSpace.inertia_data.Fmatrix * MappedInternalQuadrature.Points).colwise() + localSpace.inertia_data.translation;
    localSpace.InternalQuadrature.Weights = localSpace.inertia_data.absDetFmatrix * MappedInternalQuadrature.Weights;

#ifdef TEST
    if (abs(localSpace.InternalQuadrature.Weights.sum() - polygon.Measure) >= 1.0e-12)
        throw runtime_error("Weights inertia are wrong - 2");
#endif

    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData MappedBoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 localSpace.inertia_polygon.Vertices,
                                                 localSpace.inertia_polygon.EdgesLength,
                                                 localSpace.inertia_polygon.EdgesDirection,
                                                 localSpace.inertia_polygon.EdgesTangent,
                                                 localSpace.inertia_polygon.EdgesNormal);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 polygon.Vertices,
                                                 polygon.EdgesLength,
                                                 polygon.EdgesDirection,
                                                 polygon.EdgesTangent,
                                                 polygon.EdgesNormal);

    localSpace.Diameter = localSpace.inertia_polygon.Diameter;
    localSpace.Centroid = localSpace.inertia_polygon.Centroid;
    localSpace.Measure = localSpace.inertia_polygon.Measure;

    localSpace.constantStiff = 1.0;
    localSpace.constantMass = polygon.Measure;

    InitializeProjectorsComputation(reference_element_data,
                                    localSpace.inertia_polygon.Vertices,
                                    localSpace.inertia_polygon.Centroid,
                                    localSpace.inertia_polygon.Diameter,
                                    MappedInternalQuadrature.Points,
                                    MappedInternalQuadrature.Weights,
                                    MappedBoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    ComputePiNabla(reference_element_data,
                   localSpace.inertia_polygon.Measure,
                   localSpace.inertia_polygon.Diameter,
                   MappedInternalQuadrature.Weights,
                   MappedBoundaryQuadrature.Quadrature.Weights,
                   MappedBoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(localSpace.inertia_polygon.Measure, localSpace);

    return localSpace;
}
// ***************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputeGeometryProperties(const Gedim::GeometryUtilities &geometryUtilities,
                                                              const Utilities::Inertia_Utilities::Inertia_Data &inertia_data,
                                                              const PCC::VEM_PCC_2D_Polygon_Geometry &polygon,
                                                              PCC::VEM_PCC_2D_Polygon_Geometry &inertia_geometric_data) const
{

    inertia_geometric_data.Vertices = inertia_data.FmatrixInv * (polygon.Vertices.colwise() - inertia_data.translation);

    // Extract original cell2D geometric information
    const unsigned int &cell2DTriangulationSize = polygon.TriangulationVertices.size();
    vector<Eigen::Matrix3d> convexCell2DTriangulationPoints(cell2DTriangulationSize);
    for (unsigned int i = 0; i < cell2DTriangulationSize; i++)
    {
        if (inertia_data.signDetQ < 0)
            for (unsigned int v = 0; v < 3; v++)
                convexCell2DTriangulationPoints[i].col(v) =
                    inertia_data.FmatrixInv * (polygon.TriangulationVertices[i].col(2 - v) - inertia_data.translation);

        else
            convexCell2DTriangulationPoints[i] =
                inertia_data.FmatrixInv * (polygon.TriangulationVertices[i].colwise() - inertia_data.translation);
    }

    // compute original cell2D area and centroids
    Eigen::VectorXd convexCell2DTriangulationAreas(cell2DTriangulationSize);
    Eigen::MatrixXd convexCell2DTriangulationCentroids(3, cell2DTriangulationSize);
    for (unsigned int cct = 0; cct < cell2DTriangulationSize; cct++)
    {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) =
            geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cct], convexCell2DTriangulationAreas[cct]);
    }

    inertia_geometric_data.Measure = convexCell2DTriangulationAreas.sum();
    inertia_geometric_data.Centroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                                        convexCell2DTriangulationAreas,
                                                                        inertia_geometric_data.Measure);

    // Compute cell2D triangulation from original cell2Ds
    inertia_geometric_data.TriangulationVertices.resize(cell2DTriangulationSize);
    for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints.size(); cct++)
        inertia_geometric_data.TriangulationVertices[cct] = convexCell2DTriangulationPoints[cct];

    inertia_geometric_data.Diameter = geometryUtilities.PolygonDiameter(inertia_geometric_data.Vertices);

    inertia_geometric_data.EdgesLength = geometryUtilities.PolygonEdgeLengths(inertia_geometric_data.Vertices);
    inertia_geometric_data.EdgesTangent = geometryUtilities.PolygonEdgeTangents(inertia_geometric_data.Vertices);

    inertia_geometric_data.EdgesNormal =
        inertia_data.signDetQ * geometryUtilities.PolygonEdgeNormals(inertia_geometric_data.Vertices);

    inertia_geometric_data.EdgesDirection = polygon.EdgesDirection;
}
//****************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                    const Eigen::MatrixXd &polygonVertices,
                                                                    const Eigen::Vector3d &polygonCentroid,
                                                                    const double &polygonDiameter,
                                                                    const Eigen::MatrixXd &internalQuadraturePoints,
                                                                    const Eigen::VectorXd &internalQuadratureWeights,
                                                                    const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                                    VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumVertexBasisFunctions = numVertices;
    localSpace.NumEdgeBasisFunctions = reference_element_data.NumDofs1D * numEdges;
    localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions;
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.NumBasisFunctions =
        localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions + localSpace.NumInternalBasisFunctions;

    localSpace.NumProjectorBasisFunctions = reference_element_data.Monomials.NumMonomials;

    localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - reference_element_data.Order - 1;

    // Compute Vandermonde matrices.
    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polygonCentroid, polygonDiameter);
    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polygonDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polygonDiameter);

    // Compute mass matrix of monomials.
    localSpace.Hmatrix = localSpace.VanderInternal.transpose() * internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;

    localSpace.Qmatrix = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions, localSpace.NumProjectorBasisFunctions);
    localSpace.QmatrixInv = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions, localSpace.NumProjectorBasisFunctions);

    // Compute LLT factorization of order-1 monomials.
    // localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();
}
//****************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                   const double &polygonMeasure,
                                                   const double &polygonDiameter,
                                                   const Eigen::VectorXd &internalQuadratureWeights,
                                                   const Eigen::VectorXd &boundaryQuadratureWeights,
                                                   const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                   VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = localSpace.VanderInternalDerivatives[0].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[0] +
                         localSpace.VanderInternalDerivatives[1].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[1];
    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    localSpace.Bmatrix.setZero(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);
    // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
    localSpace.Bmatrix.leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundaryDerivatives[0].transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
        localSpace.VanderBoundaryDerivatives[1].transpose() * boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (localSpace.Order == 1)
    {
        // B_{0j} = \int_{\partial E} \phi_j
        localSpace.Bmatrix.row(0) = boundaryQuadratureWeights;
        // G_{0j} = \int_{\partial E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderBoundary.transpose() * boundaryQuadratureWeights;
    }
    else
    {
        // G_{0j} = \int_{E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderInternal.transpose() * internalQuadratureWeights;
        // Second block of B: - \int_E \Delta m_i \phi_j
        localSpace.Bmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
            (-polygonMeasure / (polygonDiameter * polygonDiameter)) *
            (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions));
        // B_{0j} = \int_{E} \phi_j (only the first internal basis
        // function has a non-zero integral)
        localSpace.Bmatrix(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) = polygonMeasure;
    }

    localSpace.PiNabla = localSpace.Gmatrix.partialPivLu().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputePolynomialsDofs(const double &polygonMeasure, VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix = MatrixXd::Zero(localSpace.NumBasisFunctions, localSpace.NumProjectorBasisFunctions);
    // boundary degrees of freedom of monomials (values at points on
    // the boundary)
    localSpace.Dmatrix.topRows(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) = localSpace.VanderBoundary;

    if (localSpace.Order > 1)
    {
        // internal degrees of freedom of monomials (scaled moments)
        localSpace.Dmatrix.bottomRows(localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polygonMeasure;
    }
}
//****************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const double &polygonMeasure,
                                                                     const double &polygonDiameter,
                                                                     const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                                     VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    localSpace.Ematrix.resize(2, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));
    localSpace.Ematrix[0].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal();

    localSpace.Ematrix[1].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() * boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (localSpace.Order > 1)
    {
        localSpace.Ematrix[0].rightCols(localSpace.NumInternalBasisFunctions) =
            -(polygonMeasure / polygonDiameter) *
            monomials.D_x(reference_element_data.Monomials).topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions);
        localSpace.Ematrix[1].rightCols(localSpace.NumInternalBasisFunctions) =
            -(polygonMeasure / polygonDiameter) *
            monomials.D_y(reference_element_data.Monomials).topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions);
    }

    localSpace.Pi0km1Der.resize(2);
    const Eigen::LLT<Eigen::MatrixXd> H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();
    localSpace.Pi0km1Der[0] = H_km1_LLT.solve(localSpace.Ematrix[0]);
    localSpace.Pi0km1Der[1] = H_km1_LLT.solve(localSpace.Ematrix[1]);
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
