#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"

#include "LAPACK_utilities.hpp"

using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace MCC
{
//****************************************************************************
VEM_MCC_Velocity_LocalSpace_Data VEM_MCC_2D_Ortho_Velocity_LocalSpace::CreateLocalSpace(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const VEM_MCC_2D_Polygon_Geometry &polygon) const
{
    VEM_MCC_Velocity_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature = quadrature.PolygonEdgesQuadrature(reference_element_data.Quadrature.ReferenceSegmentQuadrature,
                                          polygon.Vertices,
                                          polygon.EdgesLength,
                                          polygon.EdgesDirection,
                                          polygon.EdgesTangent,
                                          polygon.EdgesNormal);

    InitializeProjectorsComputation(reference_element_data,
                                    polygon.EdgesLength.size(),
                                    polygon.Centroid,
                                    polygon.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    Eigen::MatrixXd W2;
    Eigen::MatrixXd B2Nabla;
    ComputeValuesOnBoundary(polygon.Vertices,
                            polygon.EdgesNormal,
                            polygon.EdgesDirection,
                            localSpace.BoundaryQuadrature.Quadrature.Weights,
                            W2,
                            B2Nabla,
                            localSpace);

    ComputeDivergenceCoefficients(polygon.Measure, W2, localSpace);

    ComputeL2Projectors(polygon.Measure, localSpace.InternalQuadrature.Weights, B2Nabla, localSpace);

    ComputePolynomialBasisDofs(polygon.Measure, localSpace);

    ComputeStabilizationMatrix(polygon.Measure, localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_MCC_2D_Ortho_Velocity_LocalSpace::InitializeProjectorsComputation(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const unsigned int &numEdges,
    const Eigen::Vector3d &polygonCentroid,
    const double &polygonDiameter,
    const Eigen::MatrixXd &internalQuadraturePoints,
    const Eigen::VectorXd &internalQuadratureWeights,
    const Eigen::MatrixXd &boundaryQuadraturePoints,
    VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Order = reference_element_data.Order;
    localSpace.Dimension = reference_element_data.Dimension;

    localSpace.NumBoundaryBasisFunctions = reference_element_data.NumDofs1D * numEdges;

    localSpace.Nk = (reference_element_data.Order + 1) * (reference_element_data.Order + 2) / 2;
    localSpace.NkNabla = localSpace.Nk + (reference_element_data.Order + 1);

    localSpace.NumNablaInternalBasisFunctions = localSpace.Nk - 1;
    localSpace.NumBigOPlusInternalBasisFunctions = localSpace.Nk - (reference_element_data.Order + 1);
    localSpace.NumInternalBasisFunctions = localSpace.NumNablaInternalBasisFunctions
                                           + localSpace.NumBigOPlusInternalBasisFunctions;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions
                                   + localSpace.NumInternalBasisFunctions;

    const MatrixXd VanderInternalMonomials = monomials.Vander(reference_element_data.MonomialsKp1,
                                                              internalQuadraturePoints,
                                                              polygonCentroid,
                                                              polygonDiameter);

    // Compute change of polynomial basis matrix
    ChangeOfBasis(VanderInternalMonomials, internalQuadratureWeights, localSpace);

    // Compute Vandermonde matrices.
    localSpace.VanderInternalKp1 = VanderInternalMonomials * localSpace.QmatrixKp1.transpose();

    localSpace.VanderInternal = localSpace.VanderInternalKp1.leftCols(localSpace.Nk);

    const MatrixXd VanderBoundaryMonomials = monomials.Vander(
        reference_element_data.MonomialsKp1, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundaryKp1 = VanderBoundaryMonomials * localSpace.QmatrixKp1.transpose();

    localSpace.VanderBoundary = localSpace.VanderBoundaryKp1.leftCols(localSpace.Nk);

    localSpace.TkNabla.resize(localSpace.NkNabla, localSpace.Dimension * localSpace.Nk);

    MatrixXd dx = (1.0 / polygonDiameter) * monomials.D_x(reference_element_data.MonomialsKp1).leftCols(localSpace.Nk);
    MatrixXd dy = (1.0 / polygonDiameter) * monomials.D_y(reference_element_data.MonomialsKp1).leftCols(localSpace.Nk);

    dx = (localSpace.QmatrixKp1 * dx * localSpace.QmatrixInvKp1.topLeftCorner(localSpace.Nk, localSpace.Nk))
             .bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);
    dy = (localSpace.QmatrixKp1 * dy * localSpace.QmatrixInvKp1.topLeftCorner(localSpace.Nk, localSpace.Nk))
             .bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);

    localSpace.TkNabla << dx, dy;

    // L2(E)-orthogonalization process
    MatrixXd Q2;
    MatrixXd R2;
    LAPACK_utilities::MGS(localSpace.TkNabla.transpose(), Q2, R2);

    localSpace.TkNabla = Q2.transpose();
    localSpace.QmatrixTkNablaInv = (R2).transpose();
    LAPACK_utilities::inverseTri(localSpace.QmatrixTkNablaInv, localSpace.QmatrixTkNabla, 'L', 'N');

    MatrixXd V;
    VectorXd S;
    LAPACK_utilities::svd(localSpace.TkNabla, V, S);

    MatrixXd VanderInternal2k(localSpace.Dimension * localSpace.Nk,
                              localSpace.Dimension * internalQuadratureWeights.size());

    VanderInternal2k << localSpace.VanderInternal.transpose(),
        MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()),
        MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()),
        localSpace.VanderInternal.transpose();

    const MatrixXd GkNablaVanderInternal = localSpace.TkNabla * VanderInternal2k;

    localSpace.TkBigOPlus =
        V.transpose().rightCols(localSpace.Dimension * localSpace.Nk - localSpace.NkNabla).transpose();

    const MatrixXd GkBigOPlusVanderInternal = localSpace.TkBigOPlus * VanderInternal2k;

    localSpace.GkVanderInternal.resize(localSpace.Dimension * localSpace.Nk,
                                       localSpace.Dimension * internalQuadratureWeights.size());

    localSpace.GkVanderInternal << GkNablaVanderInternal, GkBigOPlusVanderInternal;

    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    VectorXd internalWeights2k(localSpace.Dimension * internalQuadratureWeights.size());
    internalWeights2k << sqrtInternalQuadratureWeights, sqrtInternalQuadratureWeights;

    const MatrixXd tempG = internalWeights2k.asDiagonal() * localSpace.GkVanderInternal.transpose();
    localSpace.Gmatrix = tempG.transpose() * tempG;
}
//****************************************************************************
void VEM_MCC_2D_Ortho_Velocity_LocalSpace::ChangeOfBasis(const Eigen::MatrixXd &VanderInternalMonomialsKp1,
                                                         const Eigen::VectorXd &internalQuadratureWeights,
                                                         VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
{

    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();

    MatrixXd Q1;
    MatrixXd R1;
    LAPACK_utilities::MGS(VanderInternalMonomialsKp1, Q1, R1);

    // L2(E)-re-orthogonalization process
    MatrixXd Q2;
    MatrixXd R2;
    MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * Q1;
    LAPACK_utilities::MGS(temp, Q2, R2);

    localSpace.Hmatrix = (Q2.transpose() * Q2).topLeftCorner(localSpace.Nk, localSpace.Nk);

    localSpace.QmatrixInvKp1 = (R2 * R1).transpose();
    LAPACK_utilities::inverseTri(localSpace.QmatrixInvKp1, localSpace.QmatrixKp1, 'L', 'N');
}
//****************************************************************************
void VEM_MCC_2D_Ortho_Velocity_LocalSpace::ComputeL2Projectors(const double &polygonMeasure,
                                                               const Eigen::VectorXd &internalQuadratureWeights,
                                                               const Eigen::MatrixXd &B2Nabla,
                                                               VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
{

    const MatrixXd HHashtagMatrix = localSpace.VanderInternalKp1.rightCols(localSpace.NkNabla).transpose() *
                                    internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;

    const MatrixXd B1Nabla = -HHashtagMatrix * localSpace.Vmatrix;

    const MatrixXd BNabla = localSpace.QmatrixTkNabla * (B1Nabla + B2Nabla);

    MatrixXd BBigOPlus = MatrixXd::Zero(localSpace.NumBigOPlusInternalBasisFunctions, localSpace.NumBasisFunctions);
    BBigOPlus.rightCols(localSpace.NumBigOPlusInternalBasisFunctions) =
        polygonMeasure * Eigen::MatrixXd::Identity(localSpace.NumBigOPlusInternalBasisFunctions,
                                                   localSpace.NumBigOPlusInternalBasisFunctions);

    localSpace.Bmatrix.resize(localSpace.Dimension * localSpace.Nk, localSpace.NumBasisFunctions);
    localSpace.Bmatrix << BNabla, BBigOPlus;

    localSpace.Pi0k = localSpace.Gmatrix.llt().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_MCC_2D_Ortho_Velocity_LocalSpace::ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                                                         const Eigen::MatrixXd &W2,
                                                                         VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
{
    MatrixXd W1 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);

    if (localSpace.Order > 0)
    {
        W1.block(1,
                 localSpace.NumBoundaryBasisFunctions,
                 localSpace.NumNablaInternalBasisFunctions,
                 localSpace.NumNablaInternalBasisFunctions) =
            -polytopeMeasure * localSpace.QmatrixTkNablaInv.topLeftCorner(localSpace.NumNablaInternalBasisFunctions,
                                                                          localSpace.NumNablaInternalBasisFunctions);
    }

    localSpace.Wmatrix = W1 + W2;
    localSpace.Vmatrix = localSpace.Hmatrix.llt().solve(localSpace.Wmatrix);
}
//****************************************************************************
void VEM_MCC_2D_Ortho_Velocity_LocalSpace::ComputeValuesOnBoundary(const Eigen::MatrixXd &polytopeVertices,
                                                                   const Eigen::MatrixXd &edgeNormals,
                                                                   const std::vector<bool> &edgeDirections,
                                                                   const Eigen::VectorXd &boundaryQuadratureWeights,
                                                                   MatrixXd &W2,
                                                                   MatrixXd &B2Nabla,
                                                                   VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
{

    const unsigned int numVertices = polytopeVertices.cols();
    const unsigned int numEdges = numVertices;

    std::vector<Eigen::VectorXd> edgeNormalsVector(localSpace.Dimension,
                                                   Eigen::VectorXd::Zero(localSpace.NumBoundaryBasisFunctions));

    Eigen::VectorXd edgeDirectionsVector(localSpace.NumBoundaryBasisFunctions);
    // offset used below to set edge-internal quadrature points and weights.
    unsigned int edgeInternalPointsOffset = 0;
    for (unsigned int i = 0; i < numEdges; ++i)
    {
        const Eigen::VectorXd &outNormalTimesAbsMapDeterminant = edgeNormals.col(i);
        const double direction = edgeDirections[i] ? 1.0 : -1.0;

        // map edge internal quadrature points
        const unsigned int numEdgeInternalQuadraturePoints = localSpace.Order + 1;

        for (unsigned int d = 0; d < localSpace.Dimension; ++d)
        {
            edgeNormalsVector[d].segment(edgeInternalPointsOffset, numEdgeInternalQuadraturePoints) =
                Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints,
                                          outNormalTimesAbsMapDeterminant[d] * direction);
        }
        edgeDirectionsVector.segment(edgeInternalPointsOffset, numEdgeInternalQuadraturePoints) =
            Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints, direction);

        edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;
    }

    W2 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);
    W2.block(0, 0, localSpace.Nk, localSpace.NumBoundaryBasisFunctions) =
        localSpace.VanderBoundary.transpose() *
        boundaryQuadratureWeights.cwiseProduct(edgeDirectionsVector).asDiagonal();

    B2Nabla = MatrixXd::Zero(localSpace.NkNabla, localSpace.NumBasisFunctions);
    B2Nabla.block(0, 0, localSpace.NkNabla, localSpace.NumBoundaryBasisFunctions) =
        localSpace.VanderBoundaryKp1.rightCols(localSpace.NkNabla).transpose() *
        boundaryQuadratureWeights.cwiseProduct(edgeDirectionsVector).asDiagonal();

    MatrixXd concatenateEdgeNormalMatrix(localSpace.Dimension * localSpace.NumBoundaryBasisFunctions,
                                         localSpace.NumBoundaryBasisFunctions);
    const MatrixXd temp1 = edgeNormalsVector[0].asDiagonal();
    const MatrixXd temp2 = edgeNormalsVector[1].asDiagonal();
    concatenateEdgeNormalMatrix << temp1, temp2;

    MatrixXd VanderBoundary2k(localSpace.Dimension * localSpace.Nk,
                              localSpace.Dimension * boundaryQuadratureWeights.size());
    VanderBoundary2k << localSpace.VanderBoundary.transpose(),
        MatrixXd::Zero(localSpace.VanderBoundary.cols(), localSpace.VanderBoundary.rows()),
        MatrixXd::Zero(localSpace.VanderBoundary.cols(), localSpace.VanderBoundary.rows()),
        localSpace.VanderBoundary.transpose();

    const MatrixXd GkNablaVanderBoundary = localSpace.TkNabla * VanderBoundary2k;

    const MatrixXd GkBigOPlusVanderBoundary = localSpace.TkBigOPlus * VanderBoundary2k;

    MatrixXd GkVanderBoundary(localSpace.Dimension * localSpace.Nk,
                              localSpace.Dimension * boundaryQuadratureWeights.size());

    GkVanderBoundary << GkNablaVanderBoundary, GkBigOPlusVanderBoundary;

    localSpace.GkVanderBoundaryTimesNormal = (GkVanderBoundary * concatenateEdgeNormalMatrix).transpose();
}
//****************************************************************************
} // namespace MCC
} // namespace VEM
} // namespace Polydim
