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

#include "VEM_PCC_3D_Inertia_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
//****************************************************************************
VEM_PCC_3D_LocalSpace_Data VEM_PCC_3D_Inertia_LocalSpace::CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                                           const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                                           const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                                           const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    VEM_PCC_3D_LocalSpace_Data localSpace;
    Quadrature::VEM_Quadrature_3D quadrature3D;

    Utilities::Inertia_Utilities inertia_utilities;
    std::vector<double> inertia_faces_measure;

    inertia_utilities.InertiaMapping3D(geometryUtilities,
                                       polyhedron.Vertices,
                                       polyhedron.Centroid,
                                       polyhedron.Diameter,
                                       polyhedron.TetrahedronVertices,
                                       localSpace.inertia_data);

    ComputeGeometryProperties(geometryUtilities, localSpace.inertia_data, polyhedron, localSpace.inertia_polyhedron, inertia_faces_measure);

    // Compute quadrature on faces
    const unsigned int numFaces = polyhedron.Faces.size();

    // Compute VEM values on faces
    VEM_PCC_2D_Inertia_LocalSpace faceVemValues;
    std::vector<Eigen::MatrixXd> facesQuadraturePoints(numFaces);
    std::vector<Eigen::VectorXd> facesQuadratureWeights(numFaces);
    std::vector<Eigen::VectorXd> facesMappedFacesQuadratureWeights(numFaces);
    std::vector<double> facesMeasure(numFaces);
    localSpace.facesLocalSpace.resize(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
    {
        localSpace.facesLocalSpace[f] = faceVemValues.Compute3DUtilities(reference_element_data_2D, polygonalFaces[f]);

        facesQuadraturePoints[f] = localSpace.facesLocalSpace[f].InternalQuadrature.Points;
        facesQuadratureWeights[f] = localSpace.facesLocalSpace[f].InternalQuadrature.Weights;
        facesMappedFacesQuadratureWeights[f] = localSpace.facesLocalSpace[f].InternalQuadrature.Weights /
                                               localSpace.facesLocalSpace[f].inertia_data.absDetFmatrix;
        facesMeasure[f] = localSpace.facesLocalSpace[f].inertia_polygon.Measure;
    }

    if (reference_element_data_2D.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
        localSpace.EdgeInternalPoints = reference_element_data_2D.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);

    localSpace.EdgeBasisCoefficients =
        utilities.ComputeEdgeBasisCoefficients(reference_element_data_2D.Order, localSpace.EdgeInternalPoints);

    Gedim::Quadrature::QuadratureData MappedInternalQuadrature =
        quadrature3D.PolyhedronInternalQuadrature(reference_element_data_3D.Quadrature,
                                                  geometryUtilities,
                                                  localSpace.inertia_polyhedron.TetrahedronVertices);

#ifdef TEST_INERTIA
    if (abs(MappedInternalQuadrature.Weights.sum() - localSpace.inertia_data.Measure) >= 1.0e-12)
        throw runtime_error("Internal Weights bulk inertia are wrong");
#endif

    localSpace.InternalQuadrature.Points =
        (localSpace.inertia_data.Fmatrix * MappedInternalQuadrature.Points).colwise() + localSpace.inertia_data.translation;
    localSpace.InternalQuadrature.Weights = localSpace.inertia_data.absDetFmatrix * MappedInternalQuadrature.Weights;

    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC MappedBoundaryQuadrature =
        quadrature3D.PolyhedronFacesQuadrature(geometryUtilities,
                                               localSpace.inertia_polyhedron.Faces,
                                               polyhedron.FacesRotationMatrix,
                                               polyhedron.FacesTranslation,
                                               localSpace.inertia_polyhedron.FacesNormal,
                                               localSpace.inertia_polyhedron.FacesNormalDirection,
                                               facesQuadraturePoints,
                                               facesMappedFacesQuadratureWeights);

    MappedBoundaryQuadrature.Quadrature.Points =
        localSpace.inertia_data.FmatrixInv *
        (MappedBoundaryQuadrature.Quadrature.Points.colwise() - localSpace.inertia_data.translation);

    unsigned int boundaryPointsOffset = 0;
    for (unsigned int f = 0; f < numFaces; f++)
    {
        const double ff = polygonalFaces[f].Measure / inertia_faces_measure[f];
        const unsigned int numFacesBoundaryPoints = facesMappedFacesQuadratureWeights[f].size();

        for (unsigned int d = 0; d < 3; d++)
            MappedBoundaryQuadrature.WeightsTimesNormal[d].segment(boundaryPointsOffset, numFacesBoundaryPoints) =
                MappedBoundaryQuadrature.WeightsTimesNormal[d].segment(boundaryPointsOffset, numFacesBoundaryPoints) *
                localSpace.facesLocalSpace[f].inertia_data.absDetFmatrix / ff;

        boundaryPointsOffset += numFacesBoundaryPoints;
    }

#ifdef TEST_INERTIA
    if (abs(localSpace.InternalQuadrature.Weights.sum() - polyhedron.Measure) >= 1.0e-12)
        throw runtime_error("Weights inertia are wrong - 2");
#endif

    localSpace.BoundaryQuadrature = quadrature3D.PolyhedronFacesQuadrature(geometryUtilities,
                                                                           polyhedron.Faces,
                                                                           polyhedron.FacesRotationMatrix,
                                                                           polyhedron.FacesTranslation,
                                                                           polyhedron.FacesNormal,
                                                                           polyhedron.FacesNormalDirection,
                                                                           facesQuadraturePoints,
                                                                           facesQuadratureWeights);

    localSpace.BoundaryQuadrature.Quadrature.Weights = MappedBoundaryQuadrature.Quadrature.Weights;

    const Eigen::MatrixXd mappedEdgeInternalQuadraturePoints =
        quadrature3D.PolyhedronInternalEdgesQuadraturePoints(reference_element_data_2D.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                             localSpace.inertia_polyhedron.Vertices,
                                                             localSpace.inertia_polyhedron.Edges,
                                                             localSpace.inertia_polyhedron.EdgesDirection,
                                                             localSpace.inertia_polyhedron.EdgesTangent);

    const MatrixXd edgeInternalQuadraturePoints =
        (localSpace.inertia_data.Fmatrix * mappedEdgeInternalQuadraturePoints).colwise() + localSpace.inertia_data.translation;

    InitializeProjectorsComputation(reference_element_data_3D,
                                    localSpace.inertia_polyhedron.Vertices,
                                    localSpace.inertia_polyhedron.Edges,
                                    localSpace.inertia_polyhedron.Faces,
                                    localSpace.inertia_polyhedron.Centroid,
                                    localSpace.inertia_polyhedron.Diameter,
                                    MappedInternalQuadrature.Points,
                                    MappedInternalQuadrature.Weights,
                                    MappedBoundaryQuadrature.Quadrature.Points,
                                    mappedEdgeInternalQuadraturePoints,
                                    localSpace);

    ComputeFaceProjectors(faceVemValues,
                          localSpace.inertia_polyhedron.Faces,
                          facesMeasure,
                          MappedBoundaryQuadrature.Quadrature.Points,
                          MappedBoundaryQuadrature.Quadrature.Weights,
                          localSpace);

    ComputePiNabla(reference_element_data_3D,
                   localSpace.inertia_polyhedron.Measure,
                   localSpace.inertia_polyhedron.Diameter,
                   MappedInternalQuadrature.Weights,
                   MappedBoundaryQuadrature.Quadrature.Weights,
                   MappedBoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(localSpace.inertia_polyhedron.Measure, localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data_3D,
                                     localSpace.inertia_polyhedron.Measure,
                                     localSpace.inertia_polyhedron.Diameter,
                                     MappedBoundaryQuadrature.WeightsTimesNormal,
                                     localSpace);

    ComputePolynomialsDofs(localSpace.inertia_polyhedron.Measure, localSpace);

    localSpace.Diameter = localSpace.inertia_polyhedron.Diameter;
    localSpace.Centroid = localSpace.inertia_polyhedron.Centroid;
    localSpace.Measure = localSpace.inertia_polyhedron.Measure;

    localSpace.constantStiff = polyhedron.Diameter;
    localSpace.constantMass = polyhedron.Measure;

    return localSpace;
}
// ***************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::ComputeGeometryProperties(const Gedim::GeometryUtilities &geometryUtilities,
                                                              const Utilities::Inertia_Data &inertia_data,
                                                              const PCC::VEM_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                              PCC::VEM_PCC_3D_Polyhedron_Geometry &inertia_geometric_data,
                                                              std::vector<double> &inertia_faces_measure) const
{

    inertia_geometric_data.Vertices = inertia_data.FmatrixInv * (polyhedron.Vertices.colwise() - inertia_data.translation);

    const std::vector<Eigen::MatrixXd> &mappedFaces3DVertices =
        geometryUtilities.PolyhedronFaceVertices(inertia_geometric_data.Vertices, polyhedron.Faces);

    inertia_geometric_data.Edges = polyhedron.Edges;
    inertia_geometric_data.Faces = polyhedron.Faces;
    inertia_geometric_data.EdgesTangent =
        geometryUtilities.PolyhedronEdgeTangents(inertia_geometric_data.Vertices, polyhedron.Edges);

    inertia_geometric_data.EdgesDirection = polyhedron.EdgesDirection;

    inertia_geometric_data.FacesTranslation = geometryUtilities.PolyhedronFaceTranslations(mappedFaces3DVertices);
    inertia_geometric_data.FacesNormal = geometryUtilities.PolyhedronFaceNormals(mappedFaces3DVertices);
    inertia_geometric_data.FacesNormalDirection =
        geometryUtilities.PolyhedronFaceNormalDirections(mappedFaces3DVertices,
                                                         geometryUtilities.PolyhedronBarycenter(inertia_geometric_data.Vertices),
                                                         inertia_geometric_data.FacesNormal);

    inertia_geometric_data.FacesRotationMatrix =
        geometryUtilities.PolyhedronFaceRotationMatrices(mappedFaces3DVertices,
                                                         inertia_geometric_data.FacesNormal,
                                                         inertia_geometric_data.FacesTranslation);

    inertia_geometric_data.Diameter = geometryUtilities.PolyhedronDiameter(inertia_geometric_data.Vertices);

    const vector<vector<unsigned int>> polyhedronFaceTriangulations =
        geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces, mappedFaces3DVertices);

    std::vector<MatrixXd> faces2DVertices =
        geometryUtilities.PolyhedronFaceRotatedVertices(mappedFaces3DVertices,
                                                        inertia_geometric_data.FacesTranslation,
                                                        inertia_geometric_data.FacesRotationMatrix);

    std::vector<std::vector<Matrix3d>> faces2DTriangulations =
        geometryUtilities.PolyhedronFaceExtractTriangulationPoints(faces2DVertices, polyhedronFaceTriangulations);

    const unsigned int numFaces = faces2DVertices.size();
    inertia_faces_measure.resize(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
    {
        // Extract original cell2D geometric information
        const vector<Eigen::Matrix3d> &convexCell2DTriangulationPoints = faces2DTriangulations[f];
        const unsigned int &numConvexCell2DTriangulation = convexCell2DTriangulationPoints.size();

        // compute original cell2D area and centroids
        Eigen::VectorXd convexCell2DTriangulationAreas(numConvexCell2DTriangulation);
        for (unsigned int cct = 0; cct < numConvexCell2DTriangulation; cct++)
            convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);

        inertia_faces_measure[f] = convexCell2DTriangulationAreas.sum();
    }

    inertia_geometric_data.Measure =
        geometryUtilities.PolyhedronVolumeByBoundaryIntegral(faces2DTriangulations,
                                                             inertia_geometric_data.FacesNormal,
                                                             inertia_geometric_data.FacesNormalDirection,
                                                             inertia_geometric_data.FacesTranslation,
                                                             inertia_geometric_data.FacesRotationMatrix);

    inertia_geometric_data.Centroid = geometryUtilities.PolyhedronCentroid(faces2DTriangulations,
                                                                           inertia_geometric_data.FacesNormal,
                                                                           inertia_geometric_data.FacesNormalDirection,
                                                                           inertia_geometric_data.FacesTranslation,
                                                                           inertia_geometric_data.FacesRotationMatrix,
                                                                           inertia_geometric_data.Measure);

    const vector<unsigned int> polyhedronTetrahedrons =
        geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(inertia_geometric_data.Vertices,
                                                                     polyhedron.Faces,
                                                                     polyhedronFaceTriangulations,
                                                                     inertia_geometric_data.Centroid);

    inertia_geometric_data.TetrahedronVertices =
        geometryUtilities.ExtractTetrahedronPoints(inertia_geometric_data.Vertices, inertia_geometric_data.Centroid, polyhedronTetrahedrons);
}
//****************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                    const Eigen::MatrixXd &polyhedronVertices,
                                                                    const Eigen::MatrixXi &polyhedronEdges,
                                                                    const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                                                    const Eigen::Vector3d &polyhedronCentroid,
                                                                    const double &polyhedronDiameter,
                                                                    const Eigen::MatrixXd &internalQuadraturePoints,
                                                                    const Eigen::VectorXd &internalQuadratureWeights,
                                                                    const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                                    const Eigen::MatrixXd &edgeInternalQuadraturePoints,
                                                                    VEM_PCC_3D_LocalSpace_Data &localSpace) const
{
    const unsigned int numVertices = polyhedronVertices.cols();
    const unsigned int numEdges = polyhedronEdges.cols();
    const unsigned int numFaces = polyhedronFaces.size();

    localSpace.Order = reference_element_data.Order;
    localSpace.Dimension = reference_element_data.Dimension;

    localSpace.NumEdgeDofs = reference_element_data.NumDofs1D;
    localSpace.NumFaceDofs = reference_element_data.NumDofs2D;
    localSpace.NumVertexBasisFunctions = numVertices;
    localSpace.NumEdgeBasisFunctions = localSpace.NumEdgeDofs * numEdges;
    localSpace.NumFaceBasisFunctions = localSpace.NumFaceDofs * numFaces;
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs3D;

    localSpace.NumBasisFunctions = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions +
                                   localSpace.NumFaceBasisFunctions + localSpace.NumInternalBasisFunctions;

    localSpace.NumProjectorBasisFunctions = (localSpace.Order + 1) * (localSpace.Order + 2) * (localSpace.Order + 3) / 6;

    localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - (localSpace.Order + 1) * (localSpace.Order + 2) / 2;

    localSpace.NumBoundaryBasisFunctions =
        localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions + localSpace.NumFaceBasisFunctions;

    // Compute Vandermonde matrices.
    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polyhedronDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polyhedronDiameter);

    // Compute positions of degrees of freedom corresponding to pointwise evaluations.
    localSpace.PointEdgeDofsCoordinates.resize(3, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions);

    localSpace.PointEdgeDofsCoordinates << polyhedronVertices, edgeInternalQuadraturePoints;

    // edge degrees of freedom of monomials (values at points on the edges).
    localSpace.VanderEdgeDofs =
        monomials.Vander(reference_element_data.Monomials, localSpace.PointEdgeDofsCoordinates, polyhedronCentroid, polyhedronDiameter);

    // Compute mass matrix of monomials.
    const VectorXd internalQuadratureWeightsSqrt = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = internalQuadratureWeightsSqrt.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;

    localSpace.Qmatrix = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions, localSpace.NumProjectorBasisFunctions);
    localSpace.QmatrixInv = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions, localSpace.NumProjectorBasisFunctions);
}
//****************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::ComputeFaceProjectors(const VEM_PCC_2D_Inertia_LocalSpace &faceVemValues,
                                                          const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                                          const std::vector<double> &facesMeasure,
                                                          const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                          const Eigen::VectorXd &boundaryQuadratureWeights,
                                                          VEM_PCC_3D_LocalSpace_Data &localSpace) const
{
    const unsigned int numFaces = polyhedronFaces.size();
    unsigned int numFaceInternalQuadrature = boundaryQuadraturePoints.cols();

    localSpace.VanderFaceProjections.setZero(numFaceInternalQuadrature, localSpace.NumBoundaryBasisFunctions);

    unsigned int faceQuadraturePointsOffset = 0;
    unsigned int faceDofsOffset = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions;

    localSpace.FaceScaledMomentsBasis.resize(numFaces);
    localSpace.FaceProjectedBasisFunctionsValues.resize(numFaces);
    for (unsigned int numFace = 0; numFace < numFaces; numFace++)
    {
        // Compute values of 2D tangential monomials at rotated quadrature points.
        const MatrixXd facePolynomialBasisValues = faceVemValues.ComputePolynomialsValues(localSpace.facesLocalSpace[numFace]);

        // Compute values of the 2D Pi^0_{k-1} projection of basis functions at rotated quadrature
        // points.
        localSpace.FaceProjectedBasisFunctionsValues[numFace] =
            faceVemValues.ComputeBasisFunctionsValues(localSpace.facesLocalSpace[numFace], ProjectionTypes::Pi0km1);

        // Fill columns of vanderFaceProjections relative to vertex basis functions.
        unsigned int numQuadraturePointsOnFace = localSpace.FaceProjectedBasisFunctionsValues[numFace].rows();
        for (unsigned int numFaceVertex = 0; numFaceVertex < polyhedronFaces[numFace].cols(); numFaceVertex++)
        {
            localSpace.VanderFaceProjections.col(polyhedronFaces[numFace](0, numFaceVertex)).segment(faceQuadraturePointsOffset, numQuadraturePointsOnFace) =
                localSpace.FaceProjectedBasisFunctionsValues[numFace].col(numFaceVertex);
        }

        if (localSpace.Order > 1)
        {
            // Compute values of tangential monomials of order k-2 divided by the area of the face. This
            // is used to compute the moments of 3D monomials on each face (see method
            // ComputePolynomialDofs()).
            localSpace.FaceScaledMomentsBasis[numFace] =
                (1.0 / facesMeasure[numFace]) *
                facePolynomialBasisValues.leftCols(localSpace.facesLocalSpace[numFace].NumInternalBasisFunctions);

            // fill columns of vanderFaceProjections relative to edge basis functions.
            for (unsigned int numFaceEdge = 0; numFaceEdge < polyhedronFaces[numFace].cols(); numFaceEdge++)
            {
                localSpace.VanderFaceProjections.block(faceQuadraturePointsOffset,
                                                       // number below is the index of the first column relative to
                                                       // this edge (assuming all edges have the same number of
                                                       // dofs)
                                                       localSpace.NumVertexBasisFunctions +
                                                           localSpace.NumEdgeDofs * polyhedronFaces[numFace](1, numFaceEdge),
                                                       numQuadraturePointsOnFace,
                                                       localSpace.NumEdgeDofs) =
                    localSpace.FaceProjectedBasisFunctionsValues[numFace].block(
                        0,
                        localSpace.facesLocalSpace[numFace].NumVertexBasisFunctions + localSpace.NumEdgeDofs * numFaceEdge,
                        numQuadraturePointsOnFace,
                        localSpace.NumEdgeDofs);
            }

            // fill columns of vanderFaceProjections relative to face internal basis functions.
            localSpace.VanderFaceProjections.block(faceQuadraturePointsOffset,
                                                   faceDofsOffset,
                                                   numQuadraturePointsOnFace,
                                                   localSpace.NumFaceDofs) =
                localSpace.FaceProjectedBasisFunctionsValues[numFace].rightCols(localSpace.NumFaceDofs);
        }

        faceQuadraturePointsOffset += numQuadraturePointsOnFace;
        faceDofsOffset += localSpace.NumFaceDofs;
    }

    localSpace.ScaledHmatrixOnBoundary.resize(localSpace.NumFaceBasisFunctions, localSpace.NumProjectorBasisFunctions);

    faceQuadraturePointsOffset = 0;
    faceDofsOffset = 0;

    if (localSpace.Order > 1)
    {
        // internal scaled moments on each face.
        for (unsigned int i = 0; i < localSpace.FaceScaledMomentsBasis.size(); ++i)
        {
            localSpace.ScaledHmatrixOnBoundary.block(faceDofsOffset,
                                                     0,
                                                     localSpace.FaceScaledMomentsBasis[i].cols(),
                                                     localSpace.NumProjectorBasisFunctions) =
                localSpace.FaceScaledMomentsBasis[i].transpose() *
                boundaryQuadratureWeights
                    .segment(faceQuadraturePointsOffset, localSpace.FaceScaledMomentsBasis[i].rows())
                    .asDiagonal() *
                localSpace.VanderBoundary.block(faceQuadraturePointsOffset,
                                                0,
                                                localSpace.FaceScaledMomentsBasis[i].rows(),
                                                localSpace.NumProjectorBasisFunctions);

            faceQuadraturePointsOffset += localSpace.FaceScaledMomentsBasis[i].rows();
            faceDofsOffset += localSpace.FaceScaledMomentsBasis[i].cols();
        }
    }
}
//****************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::ComputePiNabla(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                   const double &polyhedronMeasure,
                                                   const double &polyhedronDiameter,
                                                   const Eigen::VectorXd &internalQuadratureWeights,
                                                   const Eigen::VectorXd &boundaryQuadratureWeights,
                                                   const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                   VEM_PCC_3D_LocalSpace_Data &localSpace) const
{
    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions, localSpace.NumProjectorBasisFunctions);

    const VectorXd internalQuadratureWeightsSqrt = internalQuadratureWeights.array().sqrt();
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        const MatrixXd temp = internalQuadratureWeightsSqrt.asDiagonal() * localSpace.VanderInternalDerivatives[d];

        localSpace.Gmatrix += temp.transpose() * temp;
    }

    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    localSpace.Bmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);

    // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
    localSpace.Bmatrix.leftCols(localSpace.NumBoundaryBasisFunctions) =
        (localSpace.VanderBoundaryDerivatives[0].transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
         localSpace.VanderBoundaryDerivatives[1].transpose() * boundaryQuadratureWeightsTimesNormal[1].asDiagonal() +
         localSpace.VanderBoundaryDerivatives[2].transpose() * boundaryQuadratureWeightsTimesNormal[2].asDiagonal()) *
        localSpace.VanderFaceProjections;

    if (localSpace.Order == 1)
    {
        // B_{0j} = \int_{\partial E} \phi_j
        localSpace.Bmatrix.row(0) = localSpace.VanderFaceProjections.transpose() * boundaryQuadratureWeights;
        // G_{0j} = \int_{\partial E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderBoundary.transpose() * boundaryQuadratureWeights;
    }
    else
    {
        // G_{0j} = \int_{E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderInternal.transpose() * internalQuadratureWeights;
        // Second block of B: - \int_E \Delta m_i \phi_j
        localSpace.Bmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
            (-polyhedronMeasure / (polyhedronDiameter * polyhedronDiameter)) *
            (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions));
        // B_{0j} = \int_{E} \phi_j (only the first internal basis
        // function has a non-zero integral)
        localSpace.Bmatrix(0, localSpace.NumBasisFunctions - localSpace.NumInternalBasisFunctions) = polyhedronMeasure;
    }

    localSpace.PiNabla = localSpace.Gmatrix.partialPivLu().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::ComputePolynomialsDofs(const double &polyhedronMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix.setZero(localSpace.NumBasisFunctions, localSpace.NumProjectorBasisFunctions);

    localSpace.Dmatrix.topRows(localSpace.VanderEdgeDofs.rows()) = localSpace.VanderEdgeDofs;

    if (localSpace.Order > 1)
    {
        // internal scaled moments on each face.
        localSpace.Dmatrix.block(localSpace.VanderEdgeDofs.rows(), 0, localSpace.NumFaceBasisFunctions, localSpace.NumProjectorBasisFunctions) =
            localSpace.ScaledHmatrixOnBoundary;

        // internal scaled moments
        localSpace.Dmatrix.bottomRows(localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polyhedronMeasure;
    }
}
//****************************************************************************
void VEM_PCC_3D_Inertia_LocalSpace::ComputeL2ProjectorsOfDerivatives(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                     const double &polyhedronMeasure,
                                                                     const double &polyhedronDiameter,
                                                                     const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                                     VEM_PCC_3D_LocalSpace_Data &localSpace) const
{
    localSpace.Pi0km1Der.resize(localSpace.Dimension);

    localSpace.Ematrix.resize(3, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));
    const Eigen::LLT<Eigen::MatrixXd> H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        localSpace.Ematrix[d].leftCols(localSpace.NumBasisFunctions - localSpace.NumInternalBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
            boundaryQuadratureWeightsTimesNormal[d].asDiagonal() * localSpace.VanderFaceProjections;

        if (localSpace.Order > 1)
        {
            localSpace.Ematrix[d].rightCols(localSpace.NumInternalBasisFunctions) =
                -(polyhedronMeasure / polyhedronDiameter) *
                monomials.DerivativeMatrix(reference_element_data.Monomials, d)
                    .topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions);
        }

        localSpace.Ematrix[d] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1) * localSpace.Ematrix[d];

        localSpace.Pi0km1Der[d] = H_km1_LLT.solve(localSpace.Ematrix[d]);
    }
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
