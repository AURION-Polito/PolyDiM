#include "VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace.hpp"
#include "LAPACK_utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace MCC
{
//****************************************************************************
VEM_MCC_2D_Velocity_LocalSpace_Data VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace::CreateLocalSpace(
    const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
    const VEM_MCC_2D_Polygon_Geometry &polygon) const
{
    VEM_MCC_2D_Velocity_LocalSpace_Data localSpace;

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
                                    polygon.Measure,
                                    polygon.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    vector<MatrixXd> Cmatrix;
    utilities.MonomialTraceOnEdges(localSpace.Order,
                                   polygon.Vertices,
                                   polygon.Diameter,
                                   polygon.Centroid,
                                   polygon.EdgesDirection,
                                   polygon.EdgesTangent,
                                   Cmatrix);

    Eigen::MatrixXd W2;
    Eigen::MatrixXd B2Nabla;
    ComputeValuesOnBoundary(reference_element_data,
                            polygon.Vertices,
                            polygon.EdgesLength,
                            polygon.EdgesNormal,
                            polygon.EdgesDirection,
                            Cmatrix,
                            W2,
                            B2Nabla,
                            localSpace);

    ComputeDivergenceCoefficients(polygon.Measure, W2, localSpace);

    ComputeL2Projectors(polygon.Measure, localSpace.InternalQuadrature.Weights, B2Nabla, localSpace);

    ComputePolynomialBasisDofs(polygon.Measure, localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace::InitializeProjectorsComputation(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                               const unsigned int &numEdges,
                                                                               const Eigen::Vector3d &polygonCentroid,
                                                                               const double &polygonMeasure,
                                                                               const double &polygonDiameter,
                                                                               const Eigen::MatrixXd &internalQuadraturePoints,
                                                                               const Eigen::VectorXd &internalQuadratureWeights,
                                                                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                                               VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.Nk = (localSpace.Order + 1) * (localSpace.Order + 2) / 2;
    localSpace.NkNabla = localSpace.Nk + (localSpace.Order + 1);

    localSpace.NumBoundaryBasisFunctions = reference_element_data.NumDofs1D * numEdges;
    localSpace.NumNablaInternalBasisFunctions = localSpace.Nk - 1;
    localSpace.NumBigOPlusInternalBasisFunctions = localSpace.Nk - (localSpace.Order + 1);
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions + localSpace.NumInternalBasisFunctions;

    // Compute Vandermonde matrices
    localSpace.Centroid = polygonCentroid;
    localSpace.Diameter = polygonDiameter;
    localSpace.Measure = polygonMeasure;

    localSpace.VanderInternalKp1 =
        monomials.Vander(reference_element_data.MonomialsKp1, internalQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderInternal = localSpace.VanderInternalKp1.leftCols(localSpace.Nk);

    localSpace.VanderBoundaryKp1 =
        monomials.Vander(reference_element_data.MonomialsKp1, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundary = localSpace.VanderBoundaryKp1.leftCols(localSpace.Nk);

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;

    localSpace.QmatrixKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                               reference_element_data.MonomialsKp1.NumMonomials);
    localSpace.QmatrixInvKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                                  reference_element_data.MonomialsKp1.NumMonomials);

    localSpace.TkNabla.resize(localSpace.NkNabla, localSpace.Dimension * localSpace.Nk);

    const MatrixXd dx =
        (1.0 / polygonDiameter) *
        monomials.D_x(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);
    const MatrixXd dy =
        (1.0 / polygonDiameter) *
        monomials.D_y(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);

    localSpace.TkNabla << dx, dy;

    MatrixXd V;
    VectorXd S;
    LAPACK_utilities::svd(localSpace.TkNabla, V, S);

    MatrixXd VanderInternal2k(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * internalQuadratureWeights.size());
    VanderInternal2k << localSpace.VanderInternal.transpose(),
        MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()),
        MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()),
        localSpace.VanderInternal.transpose();

    const MatrixXd GkNablaVanderInternal = localSpace.TkNabla * VanderInternal2k;

    localSpace.TkBigOPlus = V.transpose().rightCols(localSpace.Dimension * localSpace.Nk - localSpace.NkNabla).transpose();

    const MatrixXd GkBigOPlusVanderInternal = localSpace.TkBigOPlus * VanderInternal2k;

    localSpace.GkVanderInternal.resize(localSpace.Dimension * localSpace.Nk,
                                       localSpace.Dimension * internalQuadratureWeights.size());

    localSpace.GkVanderInternal << GkNablaVanderInternal, GkBigOPlusVanderInternal;

    VectorXd internalWeights2k(localSpace.Dimension * internalQuadratureWeights.size());
    internalWeights2k << sqrtInternalQuadratureWeights, sqrtInternalQuadratureWeights;

    const MatrixXd tempG = internalWeights2k.asDiagonal() * localSpace.GkVanderInternal.transpose();
    localSpace.Gmatrix = tempG.transpose() * tempG;
}
//****************************************************************************
void VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace::ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                                                             const Eigen::MatrixXd &W2,
                                                                             VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    MatrixXd W1 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);

    if (localSpace.Order > 0)
    {
        W1.block(1, localSpace.NumBoundaryBasisFunctions, localSpace.NumNablaInternalBasisFunctions, localSpace.NumNablaInternalBasisFunctions) =
            -polytopeMeasure * Eigen::MatrixXd::Identity(localSpace.NumNablaInternalBasisFunctions,
                                                         localSpace.NumNablaInternalBasisFunctions);
    }

    localSpace.Wmatrix = W1 + W2;
    localSpace.Vmatrix = localSpace.Hmatrix.llt().solve(localSpace.Wmatrix);
}
//****************************************************************************
void VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace::ComputeL2Projectors(const double &polygonMeasure,
                                                                   const Eigen::VectorXd &internalQuadratureWeights,
                                                                   const Eigen::MatrixXd &B2Nabla,
                                                                   VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    const MatrixXd HHashtagMatrix = localSpace.VanderInternalKp1.rightCols(localSpace.NkNabla).transpose() *
                                    internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    const MatrixXd B1Nabla = -HHashtagMatrix * localSpace.Vmatrix;

    const MatrixXd BNabla = B1Nabla + B2Nabla;

    MatrixXd BBigOPlus = MatrixXd::Zero(localSpace.NumBigOPlusInternalBasisFunctions, localSpace.NumBasisFunctions);
    BBigOPlus.rightCols(localSpace.NumBigOPlusInternalBasisFunctions) =
        polygonMeasure * Eigen::MatrixXd::Identity(localSpace.NumBigOPlusInternalBasisFunctions,
                                                   localSpace.NumBigOPlusInternalBasisFunctions);

    localSpace.Bmatrix.resize(localSpace.Dimension * localSpace.Nk, localSpace.NumBasisFunctions);
    localSpace.Bmatrix << BNabla, BBigOPlus;

    localSpace.Pi0k = localSpace.Gmatrix.llt().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace::ComputeValuesOnBoundary(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                       const Eigen::MatrixXd &polytopeVertices,
                                                                       const VectorXd &edgeLengths,
                                                                       const Eigen::MatrixXd &edgeNormals,
                                                                       const std::vector<bool> &edgeDirections,
                                                                       const vector<Eigen::MatrixXd> &Cmatrixkp1,
                                                                       MatrixXd &W2,
                                                                       MatrixXd &B2Nabla,
                                                                       VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
{
    const unsigned int numVertices = polytopeVertices.cols();
    const unsigned int numEdges = numVertices;
    const unsigned int kp1 = localSpace.Order + 1;

    std::vector<Eigen::VectorXd> edgeNormalsVector(2, Eigen::VectorXd::Zero(localSpace.NumBoundaryBasisFunctions));
    Eigen::VectorXd edgeDirectionsVector(localSpace.NumBoundaryBasisFunctions);
    W2 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);
    B2Nabla = MatrixXd::Zero(localSpace.NkNabla, localSpace.NumBasisFunctions);
    MatrixXd EdgeMoments(localSpace.Nk, kp1 * numEdges);

    // offset used below to set edge-internal quadrature points and weights.
    unsigned int edgeInternalPointsOffset = 0;
    unsigned int offsetCols = 0;
    for (unsigned int i = 0; i < numEdges; ++i)
    {
        const double absMapDeterminant = edgeLengths(i);
        const Eigen::VectorXd &outNormalTimesAbsMapDeterminant = edgeNormals.col(i) * absMapDeterminant;
        const double direction = edgeDirections[i] ? 1.0 : -1.0;

        // map edge internal quadrature points
        const unsigned int numEdgeInternalQuadraturePoints = localSpace.Order + 1;

        for (unsigned int d = 0; d < 2; ++d)
        {
            edgeNormalsVector[d].segment(edgeInternalPointsOffset, numEdgeInternalQuadraturePoints) =
                Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints, outNormalTimesAbsMapDeterminant[d] * direction);
        }
        edgeDirectionsVector.segment(edgeInternalPointsOffset, numEdgeInternalQuadraturePoints) =
            Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints, direction);

        W2.block(0, offsetCols, localSpace.Nk, kp1) = Cmatrixkp1[i].topLeftCorner(localSpace.Nk, kp1) *
                                                      reference_element_data.edge_ortho.QmatrixInvKp1_1D.topLeftCorner(kp1, kp1);

        EdgeMoments.block(0, offsetCols, localSpace.Nk, kp1) =
            Cmatrixkp1[i].topLeftCorner(localSpace.Nk, kp1) *
            reference_element_data.edge_ortho.QmatrixInvKp1_1D.topLeftCorner(kp1, kp1) *
            reference_element_data.edge_ortho.Hmatrix1D.topLeftCorner(kp1, kp1);

        B2Nabla.block(0, offsetCols, localSpace.NkNabla, kp1) = Cmatrixkp1[i].bottomRows(localSpace.NkNabla) *
                                                                reference_element_data.edge_ortho.QmatrixInvKp1_1D *
                                                                reference_element_data.edge_ortho.Hmatrix1D.leftCols(kp1);

        edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;
        offsetCols += kp1;
    }

    W2.block(0, 0, localSpace.Nk, localSpace.NumBoundaryBasisFunctions) =
        W2.block(0, 0, localSpace.Nk, localSpace.NumBoundaryBasisFunctions) * edgeDirectionsVector.asDiagonal();

    B2Nabla.block(0, 0, localSpace.NkNabla, localSpace.NumBoundaryBasisFunctions) =
        B2Nabla.block(0, 0, localSpace.NkNabla, localSpace.NumBoundaryBasisFunctions) * edgeDirectionsVector.asDiagonal();

    MatrixXd concatenateEdgeNormalMatrix(localSpace.Dimension * localSpace.NumBoundaryBasisFunctions,
                                         localSpace.NumBoundaryBasisFunctions);
    const MatrixXd temp1 = edgeNormalsVector[0].asDiagonal();
    const MatrixXd temp2 = edgeNormalsVector[1].asDiagonal();
    concatenateEdgeNormalMatrix << temp1, temp2;

    MatrixXd EdgeMoments2k(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * localSpace.NumBoundaryBasisFunctions);

    EdgeMoments2k << EdgeMoments, MatrixXd::Zero(localSpace.Nk, localSpace.NumBoundaryBasisFunctions),
        MatrixXd::Zero(localSpace.Nk, localSpace.NumBoundaryBasisFunctions), EdgeMoments;

    const MatrixXd GkNablaMoments = localSpace.TkNabla * EdgeMoments2k;

    const MatrixXd GkBigOPlusMoments = localSpace.TkBigOPlus * EdgeMoments2k;

    MatrixXd GkMoments(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * localSpace.NumBoundaryBasisFunctions);
    GkMoments << GkNablaMoments, GkBigOPlusMoments;

    localSpace.GkVanderBoundaryTimesNormal = (GkMoments * concatenateEdgeNormalMatrix).transpose();
}
//****************************************************************************
} // namespace MCC
} // namespace VEM
} // namespace Polydim
