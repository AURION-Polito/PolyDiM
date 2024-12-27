#include "VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
//****************************************************************************
VEM_DF_PCC_2D_Velocity_LocalSpace_Data VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::CreateLocalSpace(
    const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_DF_PCC_2D_Velocity_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceSegmentInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceSegmentExtremaWeights,
                                                 polygon.Vertices,
                                                 polygon.EdgesLength,
                                                 polygon.EdgesDirection,
                                                 polygon.EdgesTangent,
                                                 polygon.EdgesNormal);

    localSpace.EdgesDOFs =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 polygon.Vertices,
                                                 polygon.EdgesLength,
                                                 polygon.EdgesDirection,
                                                 polygon.EdgesTangent,
                                                 polygon.EdgesNormal);

    InitializeProjectorsComputation(reference_element_data,
                                    polygon.Vertices,
                                    polygon.Centroid,
                                    polygon.Diameter,
                                    polygon.EdgesDirection,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace.EdgesDOFs.Quadrature.Points,
                                    reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
                                    reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                    localSpace);

    ComputePolynomialBasisDofs(localSpace.InternalQuadrature.Weights, localSpace);

    ComputeDivergenceCoefficients(polygon.Measure, localSpace.EdgesDOFs.WeightsTimesNormal, localSpace);

    ComputeCMatrixkm2(polygon.Diameter, localSpace.BoundaryQuadrature.WeightsTimesNormal, localSpace);

    ComputePiNabla(reference_element_data,
                   polygon.Measure,
                   polygon.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.EdgesDOFs.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(localSpace.InternalQuadrature.Weights, localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data, polygon.Diameter, localSpace.EdgesDOFs.WeightsTimesNormal, localSpace);

    ComputeStabilizationMatrix(localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::InitializeProjectorsComputation(
    const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const Eigen::MatrixXd &polygonVertices,
    const Eigen::Vector3d &polygonCentroid,
    const double &polygonDiameter,
    const std::vector<bool> &edgeDirections,
    const Eigen::MatrixXd &internalQuadraturePoints,
    const Eigen::VectorXd &internalQuadratureWeights,
    const Eigen::MatrixXd &boundaryQuadraturePoints,
    const Eigen::MatrixXd &boundaryDofQuadraturePoints,
    const Eigen::MatrixXd &referenceEdgeInternalPoints,
    const Eigen::MatrixXd &referenceEdgeDofInternalPoints,
    VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumVertexBasisFunctions = numVertices;
    localSpace.NumEdgeBasisFunctions = reference_element_data.NumDofs1D * numEdges;
    localSpace.NumBoundaryBasisFunctions = 2 * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions);

    localSpace.NumDivergenceInternalBasisFunctions = reference_element_data.NumDofs2D_Divergence;
    localSpace.NumBigOPlusInternalBasisFunctions = reference_element_data.NumDofs2D_BigOPlus;
    localSpace.NumInternalBasisFunctions = localSpace.NumDivergenceInternalBasisFunctions + localSpace.NumBigOPlusInternalBasisFunctions;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions + localSpace.NumInternalBasisFunctions;

    localSpace.NKp1 = (reference_element_data.Order + 3) * (reference_element_data.Order + 2) / 2;
    localSpace.Nk = (reference_element_data.Order + 1) * (reference_element_data.Order + 2) / 2;
    localSpace.Nkm1 = reference_element_data.Order * (reference_element_data.Order + 1) / 2;
    localSpace.Nkm2 = (reference_element_data.Order - 1) * reference_element_data.Order / 2;
    localSpace.Nkm3 = (reference_element_data.Order - 2) * (reference_element_data.Order - 1) / 2;

    localSpace.Centroid = polygonCentroid;
    localSpace.Diameter = polygonDiameter;

    // Compute Vandermonde matrices.
    const MatrixXd vanderInternalKp1 =
        monomials.Vander(reference_element_data.GBasis.monomials_data, internalQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderInternal = vanderInternalKp1.leftCols(localSpace.Nk);

    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polygonDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryDofQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polygonDiameter);

    localSpace.VanderBoundaryKp1 =
        monomials.Vander(reference_element_data.GBasis.monomials_data, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderGBigOPlus = g_basis.VanderGBigOPlus(reference_element_data.GBasis, localSpace.VanderInternal);
    localSpace.VectorDecompositionMatrices = reference_element_data.GBasis.VectorDecomposition;
    localSpace.VanderGBigOPluskm2.resize(localSpace.Dimension,
                                         MatrixXd::Zero(localSpace.VanderInternal.rows(), localSpace.NumBigOPlusInternalBasisFunctions));
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.VanderGBigOPluskm2[d] = localSpace.VanderGBigOPlus[d].leftCols(localSpace.Nkm3);

    // Compute vander edge projects Kp1
    if (referenceEdgeDofInternalPoints.rows() > 0)
        localSpace.ReferenceEdgeDofInternalPoints = referenceEdgeDofInternalPoints.row(0);
    else
        localSpace.ReferenceEdgeDofInternalPoints = Eigen::RowVectorXd();

    if (referenceEdgeInternalPoints.rows() > 0)
        localSpace.ReferenceEdgeInternalPoints = referenceEdgeInternalPoints.row(0);
    else
        localSpace.ReferenceEdgeInternalPoints = Eigen::RowVectorXd();

    localSpace.VanderKp1EdgeProjections.setZero(localSpace.NumVertexBasisFunctions + localSpace.Order * numEdges,
                                                localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions);

    localSpace.VanderKp1EdgeProjections.block(0, 0, localSpace.NumVertexBasisFunctions, localSpace.NumVertexBasisFunctions) =
        Eigen::MatrixXd::Identity(localSpace.NumVertexBasisFunctions, localSpace.NumVertexBasisFunctions);

    localSpace.EdgeBasisCoefficients =
        utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, localSpace.ReferenceEdgeDofInternalPoints);

    Eigen::MatrixXd edgeFunctionValues = utilities.ComputeValuesOnEdge(localSpace.ReferenceEdgeDofInternalPoints,
                                                                       reference_element_data.Order,
                                                                       localSpace.EdgeBasisCoefficients,
                                                                       localSpace.ReferenceEdgeInternalPoints.transpose());

    unsigned int edgeInternalPointsOffset = numVertices;
    for (unsigned int e = 0; e < numEdges; e++)
    {
        const unsigned int numEdgeInternalQuadraturePoints = localSpace.ReferenceEdgeInternalPoints.size();

        const unsigned int edgeOrigin = edgeDirections[e] ? e : (e + 1) % numEdges;
        const unsigned int edgeEnd = edgeDirections[e] ? (e + 1) % numEdges : e;

        localSpace.VanderKp1EdgeProjections.block(edgeInternalPointsOffset, edgeOrigin, numEdgeInternalQuadraturePoints, 1) =
            edgeFunctionValues.col(0);
        localSpace.VanderKp1EdgeProjections.block(edgeInternalPointsOffset, edgeEnd, numEdgeInternalQuadraturePoints, 1) =
            edgeFunctionValues.col(1);
        localSpace.VanderKp1EdgeProjections.block(edgeInternalPointsOffset,
                                                  localSpace.NumVertexBasisFunctions + e * (reference_element_data.Order - 1),
                                                  numEdgeInternalQuadraturePoints,
                                                  reference_element_data.Order - 1) =
            edgeFunctionValues.rightCols(reference_element_data.Order - 1);

        edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;
    }

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd tempKp1 = sqrtInternalQuadratureWeights.asDiagonal() * vanderInternalKp1;
    localSpace.HmatrixKp1 = tempKp1.transpose() * tempKp1;
    localSpace.Hmatrix = localSpace.HmatrixKp1.topLeftCorner(localSpace.Nk, localSpace.Nk);
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputeDivergenceCoefficients(const double &polygonMeasure,
                                                                              const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                                                              VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Wmatrix = MatrixXd::Zero(1, localSpace.NumBasisFunctions);

    localSpace.Wmatrix.row(0).segment(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        boundaryDofQuadratureWeightsTimesNormal[0];

    localSpace.Wmatrix.row(0).segment(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions,
                                      localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        boundaryDofQuadratureWeightsTimesNormal[1].transpose();

    localSpace.Vmatrix = (1.0 / polygonMeasure) * localSpace.Wmatrix;
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputePolynomialBasisDofs(const Eigen::VectorXd &internalQuadratureWeights,
                                                                           VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.NumBasisFunctions, localSpace.Nk));
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        // boundary degrees of freedom of monomials (values at points on the boundary)
        localSpace.Dmatrix[d].block(d * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
                                    0,
                                    localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions,
                                    localSpace.Nk) = localSpace.VanderBoundary;

        // Big o plus
        localSpace.Dmatrix[d].middleRows(localSpace.NumBoundaryBasisFunctions, localSpace.NumBigOPlusInternalBasisFunctions) =
            localSpace.VanderGBigOPluskm2[d].transpose() * internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    }
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputeStabilizationMatrix(VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    Eigen::MatrixXd staBmatrix = localSpace.Dmatrix[0] * localSpace.PiNabla[0] + localSpace.Dmatrix[1] * localSpace.PiNabla[1];

    staBmatrix.diagonal().array() -= 1;

    // staBmatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
    localSpace.StabMatrix = staBmatrix.transpose() * staBmatrix;
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputeCMatrixkm2(const double &polygonDiameter,
                                                                  const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                                                  VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Cmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions));

    // Boundary DOF
    localSpace.Cmatrix[0].middleCols(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        polygonDiameter * localSpace.VectorDecompositionMatrices[0][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
        localSpace.VanderBoundaryKp1.middleCols(1, localSpace.NKp1 - 1).transpose() *
        boundaryQuadratureWeightsTimesNormal[0].asDiagonal() * localSpace.VanderKp1EdgeProjections;

    localSpace.Cmatrix[0].middleCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions,
                                     localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        polygonDiameter * localSpace.VectorDecompositionMatrices[0][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
        localSpace.VanderBoundaryKp1.middleCols(1, localSpace.NKp1 - 1).transpose() *
        boundaryQuadratureWeightsTimesNormal[1].asDiagonal() * localSpace.VanderKp1EdgeProjections;

    localSpace.Cmatrix[1].middleCols(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        polygonDiameter * localSpace.VectorDecompositionMatrices[1][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
        localSpace.VanderBoundaryKp1.middleCols(1, localSpace.NKp1 - 1).transpose() *
        boundaryQuadratureWeightsTimesNormal[0].asDiagonal() * localSpace.VanderKp1EdgeProjections;

    localSpace.Cmatrix[1].middleCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions,
                                     localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        polygonDiameter * localSpace.VectorDecompositionMatrices[1][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
        localSpace.VanderBoundaryKp1.middleCols(1, localSpace.NKp1 - 1).transpose() *
        boundaryQuadratureWeightsTimesNormal[1].asDiagonal() * localSpace.VanderKp1EdgeProjections;

    // DOF big o plus
    localSpace.Cmatrix[0].middleCols(localSpace.NumBoundaryBasisFunctions, localSpace.NumBigOPlusInternalBasisFunctions) +=
        localSpace.VectorDecompositionMatrices[0][1].topLeftCorner(localSpace.Nk, localSpace.Nkm3);

    localSpace.Cmatrix[1].middleCols(localSpace.NumBoundaryBasisFunctions, localSpace.NumBigOPlusInternalBasisFunctions) +=
        localSpace.VectorDecompositionMatrices[1][1].topLeftCorner(localSpace.Nk, localSpace.Nkm3);

    // DOF divergenza
    localSpace.Cmatrix[0] += -polygonDiameter *
                             localSpace.VectorDecompositionMatrices[0][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
                             localSpace.HmatrixKp1.block(1, 0, localSpace.NKp1 - 1, 1) * localSpace.Vmatrix;

    localSpace.Cmatrix[1] += -polygonDiameter *
                             localSpace.VectorDecompositionMatrices[1][0].block(0, 1, localSpace.Nk, localSpace.NKp1 - 1) *
                             localSpace.HmatrixKp1.block(1, 0, localSpace.NKp1 - 1, 1) * localSpace.Vmatrix;

    localSpace.Cmatrixkm2.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nkm2, localSpace.NumBasisFunctions));
    localSpace.Cmatrixkm2[0] = localSpace.Cmatrix[0].topRows(localSpace.Nkm2);
    localSpace.Cmatrixkm2[1] = localSpace.Cmatrix[1].topRows(localSpace.Nkm2);
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights,
                                                                    VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    // DOF big o plus
    MatrixXd EnhancedMassMatrix = MatrixXd::Zero(localSpace.VanderGBigOPlus[0].cols(), localSpace.NumBasisFunctions);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        EnhancedMassMatrix += localSpace.VanderGBigOPlus[d].transpose() * internalQuadratureWeights.asDiagonal() *
                              localSpace.VanderInternal * localSpace.PiNabla[d];

    const MatrixXd tmpE = EnhancedMassMatrix.middleRows(localSpace.Nkm3, localSpace.Nkm1 - localSpace.Nkm3);

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.Cmatrix[d].bottomRows(localSpace.Nk - localSpace.Nkm2) +=
            localSpace.VectorDecompositionMatrices[d][1].block(localSpace.Nkm2,
                                                               localSpace.Nkm3,
                                                               localSpace.Nk - localSpace.Nkm2,
                                                               localSpace.Nkm1 - localSpace.Nkm3) *
            tmpE;

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
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputePiNabla(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                               const double &polygonMeasure,
                                                               const double &polygonDiameter,
                                                               const Eigen::VectorXd &internalQuadratureWeights,
                                                               const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                                               VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = localSpace.VanderInternalDerivatives[0].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[0] +
                         localSpace.VanderInternalDerivatives[1].transpose() * internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternalDerivatives[1];

    localSpace.Gmatrix.row(0) = (1.0 / polygonMeasure) * localSpace.VanderInternal.transpose() * internalQuadratureWeights;

    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    const MatrixXd BmatrixX =
        localSpace.VanderBoundaryDerivatives[0].transpose() * boundaryDofQuadratureWeightsTimesNormal[0].asDiagonal() +
        localSpace.VanderBoundaryDerivatives[1].transpose() * boundaryDofQuadratureWeightsTimesNormal[1].asDiagonal();

    localSpace.Bmatrix.resize(localSpace.Dimension, MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions));

    for (unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        localSpace.Bmatrix[d].row(0) = (1.0 / polygonMeasure) * localSpace.Cmatrixkm2[d].row(0);

        localSpace.Bmatrix[d].middleCols(d * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
                                         localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) += BmatrixX;

        localSpace.Bmatrix[d] += (-1.0 / (polygonDiameter * polygonDiameter)) *
                                 (reference_element_data.Monomials.Laplacian.topLeftCorner(localSpace.Nk, localSpace.Nkm2)) *
                                 localSpace.Cmatrixkm2[d];
    }

    const Eigen::PartialPivLU<Eigen::MatrixXd> Gmatrix_pivLu = localSpace.Gmatrix.partialPivLu();
    localSpace.PiNabla.resize(localSpace.Dimension);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        localSpace.PiNabla[d] = Gmatrix_pivLu.solve(localSpace.Bmatrix[d]);
}
//****************************************************************************
void VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace::ComputeL2ProjectorsOfDerivatives(
    const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const double &polygonDiameter,
    const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
    VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    const double invDiameter = 1.0 / polygonDiameter;
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
            localSpace.Ematrix[localSpace.Dimension * d1 + d2].middleCols(
                d1 * (localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
                localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) +=
                localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
                boundaryDofQuadratureWeightsTimesNormal[d2].asDiagonal();
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
