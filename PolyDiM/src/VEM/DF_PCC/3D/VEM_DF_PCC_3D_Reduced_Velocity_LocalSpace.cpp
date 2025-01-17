#include "VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
//****************************************************************************
VEM_DF_PCC_3D_Velocity_LocalSpace_Data VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::CreateLocalSpace(
    const PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
    const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data_3D,
    const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    VEM_DF_PCC_3D_Velocity_LocalSpace_Data localSpace;
    Quadrature::VEM_Quadrature_3D quadrature3D;

    // Compute quadrature on faces
    const unsigned int numFaces = polyhedron.Faces.size();

    // Compute VEM values on faces
    PCC::VEM_PCC_2D_LocalSpace faceVemValues;
    std::vector<Eigen::MatrixXd> facesQuadraturePoints(numFaces);
    std::vector<Eigen::VectorXd> facesQuadratureWeights(numFaces);
    std::vector<Eigen::MatrixXd> facesQuadraturePointsKL(numFaces);
    std::vector<Eigen::VectorXd> facesQuadratureWeightsKL(numFaces);
    localSpace.facesLocalSpace.resize(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
    {
        localSpace.facesLocalSpace[f] = faceVemValues.Compute3DUtilities_DF_PCC(reference_element_data_2D, polygonalFaces[f]);

        facesQuadraturePoints[f] = localSpace.facesLocalSpace[f].InternalQuadrature.Points;
        facesQuadratureWeights[f] = localSpace.facesLocalSpace[f].InternalQuadrature.Weights;
        facesQuadraturePointsKL[f] = localSpace.facesLocalSpace[f].InternalQuadratureKL.Points;
        facesQuadratureWeightsKL[f] = localSpace.facesLocalSpace[f].InternalQuadratureKL.Weights;
    }

    localSpace.InternalQuadrature =
        quadrature3D.PolyhedronInternalQuadrature(reference_element_data_3D.Quadrature, geometryUtilities, polyhedron.TetrahedronVertices);

    localSpace.BoundaryQuadrature = quadrature3D.PolyhedronFacesQuadrature(geometryUtilities,
                                                                           polyhedron.Faces,
                                                                           polyhedron.FacesRotationMatrix,
                                                                           polyhedron.FacesTranslation,
                                                                           polyhedron.FacesNormal,
                                                                           polyhedron.FacesNormalDirection,
                                                                           facesQuadraturePoints,
                                                                           facesQuadratureWeights);

    localSpace.BoundaryQuadratureKL = quadrature3D.PolyhedronFacesQuadrature(geometryUtilities,
                                                                             polyhedron.Faces,
                                                                             polyhedron.FacesRotationMatrix,
                                                                             polyhedron.FacesTranslation,
                                                                             polyhedron.FacesNormal,
                                                                             polyhedron.FacesNormalDirection,
                                                                             facesQuadraturePointsKL,
                                                                             facesQuadratureWeightsKL);

    const Eigen::MatrixXd edgeInternalQuadraturePoints =
        quadrature3D.PolyhedronInternalEdgesQuadraturePoints(reference_element_data_2D.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                             polyhedron.Vertices,
                                                             polyhedron.Edges,
                                                             polyhedron.EdgesDirection,
                                                             polyhedron.EdgesTangent);

    InitializeProjectorsComputation(reference_element_data_3D,
                                    polyhedron.Vertices,
                                    polyhedron.Edges,
                                    polyhedron.Faces,
                                    polyhedron.Centroid,
                                    polyhedron.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace.BoundaryQuadratureKL.Quadrature.Points,
                                    edgeInternalQuadraturePoints,
                                    localSpace);

    ComputeFaceProjectors(faceVemValues,
                          localSpace.facesLocalSpace,
                          polyhedron.Faces,
                          polygonalFaces,
                          polyhedron.FacesTangents,
                          polyhedron.FacesTangentsGlobalDirection,
                          polyhedron.FacesNormal,
                          polyhedron.FacesNormalDirection,
                          polyhedron.FacesNormalGlobalDirection,
                          localSpace.BoundaryQuadrature.Quadrature.Points,
                          localSpace.BoundaryQuadratureKL.Quadrature.Points,
                          localSpace);

    ComputePolynomialBasisDofs(polyhedron.Measure,
                               polyhedron.Diameter,
                               localSpace.InternalQuadrature.Weights,
                               localSpace.BoundaryQuadrature.Quadrature.Weights,
                               localSpace);

    ComputeDivergenceCoefficients(reference_element_data_3D,
                                  polyhedron.Measure,
                                  polyhedron.Diameter,
                                  polyhedron.FacesNormalGlobalDirection,
                                  polygonalFaces,
                                  localSpace);

    ComputeCMatrixkm2(polyhedron.Measure, polyhedron.Diameter, localSpace.BoundaryQuadratureKL.Quadrature.Weights, localSpace);

    ComputePiNabla(reference_element_data_3D,
                   polyhedron.Measure,
                   polyhedron.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.BoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(localSpace.InternalQuadrature.Weights, localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data_3D, polyhedron.Diameter, localSpace.BoundaryQuadrature.WeightsTimesNormal, localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::InitializeProjectorsComputation(
    const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
    const Eigen::MatrixXd &polyhedronVertices,
    const Eigen::MatrixXi &polyhedronEdges,
    const std::vector<Eigen::MatrixXi> &polyhedronFaces,
    const Eigen::Vector3d &polyhedronCentroid,
    const double &polyhedronDiameter,
    const Eigen::MatrixXd &internalQuadraturePoints,
    const Eigen::VectorXd &internalQuadratureWeights,
    const Eigen::MatrixXd &boundaryQuadraturePoints,
    const Eigen::MatrixXd &boundaryQuadratureKLPoints,
    const Eigen::MatrixXd &edgeInternalQuadraturePoints,
    VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    const unsigned int numVertices = polyhedronVertices.cols();
    const unsigned int numEdges = polyhedronEdges.cols();
    const unsigned int numFaces = polyhedronFaces.size();

    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumVertexBasisFunctions = numVertices;
    localSpace.NumEdgeBasisFunctions = reference_element_data.NumDofs1D * numEdges;
    localSpace.NumNormalBasisFunctions = reference_element_data.NumDofs2D * numFaces;
    localSpace.NumTangentsBasisFunctions = reference_element_data.NumDofs2D * 2 * numFaces;
    localSpace.NumBoundaryBasisFunctions =
        localSpace.Dimension * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) +
        localSpace.NumTangentsBasisFunctions + localSpace.NumNormalBasisFunctions;

    localSpace.NumDivergenceInternalBasisFunctions = reference_element_data.NumDofs3D_Divergence;
    localSpace.NumBigOPlusInternalBasisFunctions = reference_element_data.NumDofs3D_BigOPlus;
    localSpace.NumInternalBasisFunctions = localSpace.NumDivergenceInternalBasisFunctions + localSpace.NumBigOPlusInternalBasisFunctions;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions + localSpace.NumInternalBasisFunctions;

    localSpace.NKp1 =
        (reference_element_data.Order + 3) * (reference_element_data.Order + 2) * (reference_element_data.Order + 4) / 6;
    localSpace.Nk =
        (reference_element_data.Order + 1) * (reference_element_data.Order + 2) * (reference_element_data.Order + 3) / 6;
    localSpace.Nkm1 = reference_element_data.Order * (reference_element_data.Order + 1) * (reference_element_data.Order + 2) / 6;
    localSpace.Nkm2 = (reference_element_data.Order - 1) * reference_element_data.Order * (reference_element_data.Order + 1) / 6;
    localSpace.Nkm3 = (reference_element_data.Order - 2) * (reference_element_data.Order - 1) * reference_element_data.Order / 6;
    localSpace.Nkm4 =
        (reference_element_data.Order - 3) * (reference_element_data.Order - 2) * (reference_element_data.Order - 1) / 6;

    // Compute Vandermonde matrices.
    localSpace.Diameter = polyhedronDiameter;
    localSpace.Centroid = polyhedronCentroid;

    const MatrixXd vanderInternalKp1 =
        monomials.Vander(reference_element_data.GBasis.monomials_data, internalQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderInternal = vanderInternalKp1.leftCols(localSpace.Nk);

    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polyhedronDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polyhedronDiameter);

    localSpace.VanderBoundaryKL =
        monomials.Vander(reference_element_data.GBasis.monomials_data, boundaryQuadratureKLPoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderGBigOPlus = g_basis.VanderGBigOPlus(reference_element_data.GBasis, localSpace.VanderInternal);
    localSpace.VectorDecompositionMatrices = g_basis.VectorDecomposition(reference_element_data.GBasis);
    localSpace.VanderGBigOPluskm2.resize(3, MatrixXd::Zero(localSpace.VanderInternal.rows(), localSpace.NumBigOPlusInternalBasisFunctions));
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        localSpace.VanderGBigOPluskm2[d] << localSpace.VanderGBigOPlus[d].leftCols(localSpace.Nkm3 - localSpace.Nkm4),
            localSpace.VanderGBigOPlus[d].middleCols(localSpace.Nkm1 - localSpace.Nkm2, localSpace.Nkm3),
            localSpace.VanderGBigOPlus[d].middleCols(2 * localSpace.Nkm1 - localSpace.Nkm2, localSpace.Nkm3);
    }

    // Compute positions of degrees of freedom corresponding to pointwise evaluations.
    Eigen::MatrixXd PointEdgeDofsCoordinates =
        Eigen::MatrixXd(3, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions);

    PointEdgeDofsCoordinates << polyhedronVertices, edgeInternalQuadraturePoints;

    // edge degrees of freedom of monomials (values at points on the edges).
    localSpace.VanderEdgeDofs =
        monomials.Vander(reference_element_data.Monomials, PointEdgeDofsCoordinates, polyhedronCentroid, polyhedronDiameter);

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd tempKp1 = sqrtInternalQuadratureWeights.asDiagonal() * vanderInternalKp1;
    localSpace.HmatrixKp1 = tempKp1.transpose() * tempKp1;
    localSpace.Hmatrix = localSpace.HmatrixKp1.topLeftCorner(localSpace.Nk, localSpace.Nk);
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputeFaceProjectors(const PCC::VEM_PCC_2D_LocalSpace &faceVemValues,
                                                                      const vector<PCC::VEM_PCC_2D_LocalSpace_Data> &facesLocalSpace,
                                                                      const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                                                      const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                                      const std::vector<std::array<Eigen::Vector3d, 2>> &facesTangents,
                                                                      const std::vector<std::array<bool, 2>> &facesTangentsGlobalDirection,
                                                                      const vector<Eigen::Vector3d> &facesNormals,
                                                                      const std::vector<bool> &facesNormalDirections,
                                                                      const std::vector<bool> &facesNormalGlobalDirections,
                                                                      const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                                      const Eigen::MatrixXd &boundaryQuadratureKLPoints,
                                                                      VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    const unsigned int numFaces = polyhedronFaces.size();
    const unsigned int numFaceInternalQuadratureKL = boundaryQuadratureKLPoints.cols();
    const unsigned int numFaceInternalQuadrature = boundaryQuadraturePoints.cols();
    const unsigned int Nkm2_2D = 0.5 * localSpace.Order * (localSpace.Order - 1);

    localSpace.VanderFaceProjectionsKp1TimesNormal.setZero(numFaceInternalQuadratureKL, localSpace.NumBoundaryBasisFunctions);
    localSpace.VanderFaceProjectionsKm1.resize(localSpace.Dimension,
                                               MatrixXd::Zero(numFaceInternalQuadrature, localSpace.NumBoundaryBasisFunctions));
    localSpace.ScaledHmatrixOnBoundary.resize(
        localSpace.Dimension,
        MatrixXd::Zero(numFaceInternalQuadrature, localSpace.NumTangentsBasisFunctions + localSpace.NumNormalBasisFunctions));

    unsigned int faceQuadraturePointsOffset = 0;
    unsigned int faceQuadraturePointsOffsetKL = 0;
    const unsigned int numVerteXEdgeBasisFunctions = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions;
    unsigned int faceTangent1DofsOffset = 3 * numVerteXEdgeBasisFunctions;
    unsigned int faceTangent2DofsOffset = 3 * numVerteXEdgeBasisFunctions + Nkm2_2D * numFaces;
    unsigned int faceNormalDofsOffset = faceTangent2DofsOffset + Nkm2_2D * numFaces;
    for (unsigned int numFace = 0; numFace < numFaces; numFace++)
    {
        const double faceNormalDirection = facesNormalDirections[numFace] ? 1.0 : -1.0;
        const double faceNormalGlobalDirection = facesNormalGlobalDirections[numFace] ? 1.0 : -1.0;
        const double faceTangent1GlobalDirection = facesTangentsGlobalDirection[numFace][0] ? 1.0 : -1.0;
        const double faceTangent2GlobalDirection = facesTangentsGlobalDirection[numFace][1] ? 1.0 : -1.0;

        // Compute values of 2D tangential monomials at rotated quadrature points.
        const MatrixXd facePolynomialBasisValues =
            (1.0 / polygonalFaces[numFace].Measure) *
            faceVemValues.ComputePolynomialsValues(facesLocalSpace[numFace]).leftCols(Nkm2_2D);

        // Compute values of the 2D Pi^0_{k-1} projection of basis functions at rotated quadrature points.
        const MatrixXd thisFaceProjectedBasisFunctionsValuesKL =
            faceVemValues.ComputeBasisFunctionsValues(facesLocalSpace[numFace], PCC::ProjectionTypes::Pi0klm1);

        const MatrixXd thisFaceProjectedBasisFunctionsValues =
            faceVemValues.ComputeBasisFunctionsValues(facesLocalSpace[numFace], PCC::ProjectionTypes::Pi0km1);

        // Fill columns of vanderFaceProjections relative to vertex basis functions.
        const unsigned int numQuadraturePointsOnFaceKL = thisFaceProjectedBasisFunctionsValuesKL.rows();
        const unsigned int numQuadraturePointsOnFace = thisFaceProjectedBasisFunctionsValues.rows();
        for (unsigned int numFaceVertex = 0; numFaceVertex < polyhedronFaces[numFace].cols(); numFaceVertex++)
        {
            for (unsigned int d = 0; d < localSpace.Dimension; d++)
            {
                localSpace.VanderFaceProjectionsKp1TimesNormal
                    .col(d * numVerteXEdgeBasisFunctions + polyhedronFaces[numFace](0, numFaceVertex))
                    .segment(faceQuadraturePointsOffsetKL, numQuadraturePointsOnFaceKL) =
                    thisFaceProjectedBasisFunctionsValuesKL.col(numFaceVertex) * faceNormalDirection * facesNormals[numFace](d);

                localSpace.VanderFaceProjectionsKm1[d]
                    .col(d * numVerteXEdgeBasisFunctions + polyhedronFaces[numFace](0, numFaceVertex))
                    .segment(faceQuadraturePointsOffset, numQuadraturePointsOnFace) =
                    thisFaceProjectedBasisFunctionsValues.col(numFaceVertex);
            }
        }

        // fill columns of vanderFaceProjections relative to edge basis functions.
        for (unsigned int numFaceEdge = 0; numFaceEdge < polyhedronFaces[numFace].cols(); numFaceEdge++)
        {
            // number below is the index of the first column relative to
            // this edge (assuming all edges have the same number of
            // dofs)
            for (unsigned int d = 0; d < localSpace.Dimension; d++)
            {
                localSpace.VanderFaceProjectionsKp1TimesNormal.block(
                    faceQuadraturePointsOffsetKL,
                    d * numVerteXEdgeBasisFunctions + localSpace.NumVertexBasisFunctions +
                        (localSpace.Order - 1) * polyhedronFaces[numFace](1, numFaceEdge),
                    numQuadraturePointsOnFaceKL,
                    (localSpace.Order - 1)) =
                    thisFaceProjectedBasisFunctionsValuesKL.block(0,
                                                                  facesLocalSpace[numFace].NumVertexBasisFunctions +
                                                                      (localSpace.Order - 1) * numFaceEdge,
                                                                  numQuadraturePointsOnFaceKL,
                                                                  (localSpace.Order - 1)) *
                    faceNormalDirection * facesNormals[numFace](d);

                localSpace.VanderFaceProjectionsKm1[d].block(faceQuadraturePointsOffset,
                                                             d * numVerteXEdgeBasisFunctions + localSpace.NumVertexBasisFunctions +
                                                                 (localSpace.Order - 1) * polyhedronFaces[numFace](1, numFaceEdge),
                                                             numQuadraturePointsOnFace,
                                                             (localSpace.Order - 1)) =
                    thisFaceProjectedBasisFunctionsValues.block(0,
                                                                facesLocalSpace[numFace].NumVertexBasisFunctions +
                                                                    (localSpace.Order - 1) * numFaceEdge,
                                                                numQuadraturePointsOnFace,
                                                                (localSpace.Order - 1));
            }
        }

        // fill columns of vanderFaceProjections relative to face internal basis functions.
        localSpace.VanderFaceProjectionsKp1TimesNormal.block(faceQuadraturePointsOffsetKL, faceNormalDofsOffset, numQuadraturePointsOnFaceKL, Nkm2_2D) =
            faceNormalGlobalDirection * thisFaceProjectedBasisFunctionsValuesKL.rightCols(Nkm2_2D);

        for (unsigned int d = 0; d < localSpace.Dimension; d++)
        {
            localSpace.VanderFaceProjectionsKm1[d].block(faceQuadraturePointsOffset, faceNormalDofsOffset, numQuadraturePointsOnFace, Nkm2_2D) =
                faceNormalGlobalDirection * thisFaceProjectedBasisFunctionsValues.rightCols(Nkm2_2D) *
                faceNormalDirection * facesNormals[numFace](d);

            localSpace.VanderFaceProjectionsKm1[d].block(faceQuadraturePointsOffset, faceTangent1DofsOffset, numQuadraturePointsOnFace, Nkm2_2D) =
                thisFaceProjectedBasisFunctionsValues.rightCols(Nkm2_2D) * faceTangent1GlobalDirection *
                facesTangents[numFace][0](d);

            localSpace.VanderFaceProjectionsKm1[d].block(faceQuadraturePointsOffset, faceTangent2DofsOffset, numQuadraturePointsOnFace, Nkm2_2D) =
                thisFaceProjectedBasisFunctionsValues.rightCols(Nkm2_2D) * faceTangent2GlobalDirection *
                facesTangents[numFace][1](d);

            localSpace.ScaledHmatrixOnBoundary[d].block(faceQuadraturePointsOffset,
                                                        faceNormalDofsOffset - 3 * numVerteXEdgeBasisFunctions,
                                                        numQuadraturePointsOnFace,
                                                        Nkm2_2D) =
                faceNormalGlobalDirection * facePolynomialBasisValues * faceNormalDirection * facesNormals[numFace](d);

            localSpace.ScaledHmatrixOnBoundary[d].block(faceQuadraturePointsOffset,
                                                        faceTangent1DofsOffset - 3 * numVerteXEdgeBasisFunctions,
                                                        numQuadraturePointsOnFace,
                                                        Nkm2_2D) =
                facePolynomialBasisValues * faceTangent1GlobalDirection * facesTangents[numFace][0](d);

            localSpace.ScaledHmatrixOnBoundary[d].block(faceQuadraturePointsOffset,
                                                        faceTangent2DofsOffset - 3 * numVerteXEdgeBasisFunctions,
                                                        numQuadraturePointsOnFace,
                                                        Nkm2_2D) =
                facePolynomialBasisValues * faceTangent2GlobalDirection * facesTangents[numFace][1](d);
        }

        faceQuadraturePointsOffsetKL += numQuadraturePointsOnFaceKL;
        faceQuadraturePointsOffset += numQuadraturePointsOnFace;
        faceNormalDofsOffset += Nkm2_2D;
        faceTangent1DofsOffset += Nkm2_2D;
        faceTangent2DofsOffset += Nkm2_2D;
    }
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputeDivergenceCoefficients(
    const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
    const double &polyhedronMeasure,
    const double &polyhedronDiameter,
    const std::vector<bool> &faceNormalGlobalDirections,
    const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
    VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Wmatrix = MatrixXd::Zero(1, localSpace.NumBasisFunctions);

    const unsigned int Nkm2_2D = reference_element_data.NumDofs2D;
    unsigned int offsetDof =
        3 * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) + localSpace.NumTangentsBasisFunctions;
    for (unsigned int f = 0; f < faceNormalGlobalDirections.size(); f++)
    {
        localSpace.Wmatrix(0, offsetDof) = faceNormalGlobalDirections[f] ? polygonalFaces[f].Measure
                                                                         : -polygonalFaces[f].Measure;
        offsetDof += Nkm2_2D;
    }

    localSpace.Vmatrix = (1.0 / polyhedronMeasure) * localSpace.Wmatrix;
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputePolynomialBasisDofs(const double &polyhedronMeasure,
                                                                           const double &polyhedronDiameter,
                                                                           const Eigen::VectorXd &internalQuadratureWeights,
                                                                           const Eigen::VectorXd &boundaryQuadratureWeights,
                                                                           VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.NumBasisFunctions, localSpace.Nk));

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        // Vertex + Edge
        localSpace.Dmatrix[d].middleRows(d * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
                                         localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderEdgeDofs;

        // Face Tangent and Normal
        localSpace.Dmatrix[d].middleRows(3 * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
                                         localSpace.NumNormalBasisFunctions + localSpace.NumTangentsBasisFunctions) =
            localSpace.ScaledHmatrixOnBoundary[d].transpose() * boundaryQuadratureWeights.asDiagonal() * localSpace.VanderBoundary;

        // Big o plus
        localSpace.Dmatrix[d].middleRows(localSpace.NumBoundaryBasisFunctions, localSpace.NumBigOPlusInternalBasisFunctions) =
            (1.0 / polyhedronMeasure) * localSpace.VanderGBigOPluskm2[d].transpose() *
            internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    }
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputeCMatrixkm2(const double &polyhedronMeasure,
                                                                  const double &polyhedronDiameter,
                                                                  const Eigen::VectorXd &boundaryQuadratureKLWeights,
                                                                  VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Cmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions));

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        // Boundary DOF
        localSpace.Cmatrix[d].leftCols(localSpace.NumBoundaryBasisFunctions) =
            polyhedronDiameter * localSpace.VectorDecompositionMatrices[d][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
            localSpace.VanderBoundaryKL.middleCols(1, localSpace.NKp1 - 1).transpose() *
            boundaryQuadratureKLWeights.asDiagonal() * localSpace.VanderFaceProjectionsKp1TimesNormal;

        // DOF big o plus
        MatrixXd tmpVD = MatrixXd::Zero(localSpace.Nk, localSpace.NumBigOPlusInternalBasisFunctions);
        tmpVD << localSpace.VectorDecompositionMatrices[d][1].middleCols(0, localSpace.Nkm3 - localSpace.Nkm4),
            localSpace.VectorDecompositionMatrices[d][1].middleCols(localSpace.Nkm1 - localSpace.Nkm2, localSpace.Nkm3),
            localSpace.VectorDecompositionMatrices[d][1].middleCols(2 * localSpace.Nkm1 - localSpace.Nkm2, localSpace.Nkm3);

        localSpace.Cmatrix[d].middleCols(localSpace.NumBoundaryBasisFunctions, localSpace.NumBigOPlusInternalBasisFunctions) +=
            polyhedronMeasure * tmpVD;

        // divergence term
        localSpace.Cmatrix[d] += -polyhedronDiameter *
                                 localSpace.VectorDecompositionMatrices[d][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
                                 localSpace.HmatrixKp1.block(1, 0, localSpace.NKp1 - 1, 1) * localSpace.Vmatrix;
    }

    localSpace.Cmatrixkm2.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nkm2, localSpace.NumBasisFunctions));
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.Cmatrixkm2[d] = localSpace.Cmatrix[d].topRows(localSpace.Nkm2);
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights,
                                                                    VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    // DOF big o plus
    MatrixXd EnhancedMassMatrix = MatrixXd::Zero(localSpace.VanderGBigOPlus[0].cols(), localSpace.NumBasisFunctions);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        EnhancedMassMatrix += localSpace.VanderGBigOPlus[d].transpose() * internalQuadratureWeights.asDiagonal() *
                              localSpace.VanderInternal * localSpace.PiNabla[d];

    MatrixXd tmpE = MatrixXd::Zero(3 * (localSpace.Nkm1 - localSpace.Nkm3) + localSpace.Nkm4 - localSpace.Nkm2,
                                   localSpace.NumBasisFunctions);
    tmpE << EnhancedMassMatrix.middleRows(localSpace.Nkm3 - localSpace.Nkm4,
                                          localSpace.Nkm1 - localSpace.Nkm3 + localSpace.Nkm4 - localSpace.Nkm2),
        EnhancedMassMatrix.middleRows(localSpace.Nkm1 - localSpace.Nkm2 + localSpace.Nkm3, localSpace.Nkm1 - localSpace.Nkm3),
        EnhancedMassMatrix.middleRows(2 * localSpace.Nkm1 - localSpace.Nkm2 + localSpace.Nkm3,
                                      localSpace.Nkm1 - localSpace.Nkm3);

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        MatrixXd tmpVD = MatrixXd::Zero(localSpace.Nk - localSpace.Nkm2,
                                        3 * (localSpace.Nkm1 - localSpace.Nkm3) + localSpace.Nkm4 - localSpace.Nkm2);
        tmpVD << localSpace.VectorDecompositionMatrices[d][1].block(localSpace.Nkm2,
                                                                    localSpace.Nkm3 - localSpace.Nkm4,
                                                                    localSpace.Nk - localSpace.Nkm2,
                                                                    localSpace.Nkm1 - localSpace.Nkm3 +
                                                                        localSpace.Nkm4 - localSpace.Nkm2),
            localSpace.VectorDecompositionMatrices[d][1].block(localSpace.Nkm2,
                                                               localSpace.Nkm1 - localSpace.Nkm2 + localSpace.Nkm3,
                                                               localSpace.Nk - localSpace.Nkm2,
                                                               localSpace.Nkm1 - localSpace.Nkm3),
            localSpace.VectorDecompositionMatrices[d][1].block(localSpace.Nkm2,
                                                               2 * localSpace.Nkm1 - localSpace.Nkm2 + localSpace.Nkm3,
                                                               localSpace.Nk - localSpace.Nkm2,
                                                               localSpace.Nkm1 - localSpace.Nkm3);

        localSpace.Cmatrix[d].bottomRows(localSpace.Nk - localSpace.Nkm2) += tmpVD * tmpE;
    }

    const Eigen::LLT<Eigen::MatrixXd> Hmatrix_llt = localSpace.Hmatrix.llt();
    localSpace.Pi0k.resize(localSpace.Dimension);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.Pi0k[d] = Hmatrix_llt.solve(localSpace.Cmatrix[d]);

    const Eigen::LLT<Eigen::MatrixXd> Hmatrix_km2_llt =
        localSpace.Hmatrix.topLeftCorner(localSpace.Nkm2, localSpace.Nkm2).llt();
    localSpace.Pi0km2.resize(localSpace.Dimension);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.Pi0km2[d] = Hmatrix_km2_llt.solve(localSpace.Cmatrixkm2[d]);
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputePiNabla(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                               const double &polyhedronMeasure,
                                                               const double &polyhedronDiameter,
                                                               const Eigen::VectorXd &internalQuadratureWeights,
                                                               const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                               VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = localSpace.VanderInternalDerivatives[0].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[0] +
                         localSpace.VanderInternalDerivatives[1].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[1] +
                         localSpace.VanderInternalDerivatives[2].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[2];

    localSpace.Gmatrix.row(0) = (1.0 / polyhedronMeasure) * localSpace.VanderInternal.transpose() * internalQuadratureWeights;

    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    const MatrixXd BmatrixX =
        localSpace.VanderBoundaryDerivatives[0].transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
        localSpace.VanderBoundaryDerivatives[1].transpose() * boundaryQuadratureWeightsTimesNormal[1].asDiagonal() +
        localSpace.VanderBoundaryDerivatives[2].transpose() * boundaryQuadratureWeightsTimesNormal[2].asDiagonal();

    localSpace.Bmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions));

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        localSpace.Bmatrix[d].row(0) = (1.0 / polyhedronMeasure) * localSpace.Cmatrixkm2[d].row(0);
        localSpace.Bmatrix[d].leftCols(localSpace.NumBoundaryBasisFunctions) += BmatrixX * localSpace.VanderFaceProjectionsKm1[d];

        localSpace.Bmatrix[d] += (-1.0 / (polyhedronDiameter * polyhedronDiameter)) *
                                 (reference_element_data.Monomials.Laplacian.topLeftCorner(localSpace.Nk, localSpace.Nkm2)) *
                                 localSpace.Cmatrixkm2[d];
    }

    const Eigen::PartialPivLU<Eigen::MatrixXd> Gmatrix_pivLu = localSpace.Gmatrix.partialPivLu();
    localSpace.PiNabla.resize(localSpace.Dimension);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.PiNabla[d] = Gmatrix_pivLu.solve(localSpace.Bmatrix[d]);
}
//****************************************************************************
void VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace::ComputeL2ProjectorsOfDerivatives(
    const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
    const double &polyhedronDiameter,
    const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
    VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    const double invDiameter = 1.0 / polyhedronDiameter;
    localSpace.Ematrix.resize(localSpace.Dimension * localSpace.Dimension,
                              Eigen::MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));
    for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
    {
        for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
        {
            localSpace.Ematrix[localSpace.Dimension * d1 + d2] =
                -invDiameter *
                monomials.DerivativeMatrix(reference_element_data.Monomials, d2)
                    .topLeftCorner(localSpace.Nkm1, localSpace.Nkm2) *
                localSpace.Cmatrixkm2[d1];
            localSpace.Ematrix[localSpace.Dimension * d1 + d2].leftCols(localSpace.NumBoundaryBasisFunctions) +=
                localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
                boundaryQuadratureWeightsTimesNormal[d2].asDiagonal() * localSpace.VanderFaceProjectionsKm1[d1];
        }
    }

    localSpace.Pi0km1Der.resize(localSpace.Dimension * localSpace.Dimension);
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();

    for (unsigned int i = 0; i < localSpace.Dimension * localSpace.Dimension; i++)
        localSpace.Pi0km1Der[i] = H_km1_LLT.solve(localSpace.Ematrix[i]);
}
//****************************************************************************
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim
