#include "VEM_MCC_2D_VelocityLocalSpace.hpp"
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
VEM_MCC_VelocityLocalSpace_Data VEM_MCC_2D_VelocityLocalSpace::CreateLocalSpace(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                                const VEM_MCC_2D_Polygon_Geometry& polygon) const
{
    VEM_MCC_VelocityLocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature = quadrature.PolygonEdgesQuadrature(reference_element_data.Quadrature,
                                                                      polygon.Vertices,
                                                                      polygon.EdgesLength,
                                                                      polygon.EdgesDirection,
                                                                      polygon.EdgesTangent,
                                                                      polygon.EdgesNormal);


    InitializeProjectorsComputation(reference_element_data,
                                    polygon.Vertices,
                                    polygon.Centroid,
                                    polygon.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    ComputeL2Projectors(polygon.Measure,
                        localSpace);

    localSpace.Dmatrix = ComputePolynomialBasisDofs(polygon.Measure,
                                                    localSpace);

    localSpace.StabMatrix = ComputeStabilizationMatrixPi0k(polygon.Measure,
                                                           localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_MCC_2D_VelocityLocalSpace::InitializeProjectorsComputation(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                    const Eigen::MatrixXd& polygonVertices,
                                                                    const Eigen::Vector3d& polygonCentroid,
                                                                    const double& polygonDiameter,
                                                                    const Eigen::MatrixXd& internalQuadraturePoints,
                                                                    const Eigen::VectorXd& internalQuadratureWeights,
                                                                    const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                                    VEM_MCC_VelocityLocalSpace_Data& localSpace) const
{
    const unsigned int numVertices = polygonVertices.cols();
    const unsigned int numEdges = numVertices;

    localSpace.Order = reference_element_data.Order;

    localSpace.Nk = (localSpace.Order  + 1) * (localSpace.Order  + 2 ) / 2;
    localSpace.NkNabla = localSpace.Nk + (localSpace.Order + 1);

    localSpace.NumBoundaryBasisFunctions = reference_element_data.NumDofs1D * numEdges;
    localSpace.NumNablaInternalBasisFunctions = localSpace.Nk - 1 ;
    localSpace.NumBigOPlusInternalBasisFunctions = localSpace.Nk - (localSpace.Order + 1 );
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions +
                                   localSpace.NumInternalBasisFunctions;


    // Compute Vandermonde matrices.
    localSpace.VanderInternalKp1 = monomials.Vander(reference_element_data.MonomialsKp1,
                                                    internalQuadraturePoints,
                                                    polygonCentroid,
                                                    polygonDiameter);

    localSpace.VanderInternal = localSpace.VanderInternalKp1.leftCols(localSpace.Nk);

    localSpace.VanderBoundaryKp1 = monomials.Vander(reference_element_data.MonomialsKp1,
                                                    boundaryQuadraturePoints,
                                                    polygonCentroid,
                                                    polygonDiameter);

    localSpace.VanderBoundary = localSpace.VanderBoundaryKp1.leftCols(localSpace.Nk);


    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;

    localSpace.QmatrixKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                               reference_element_data.MonomialsKp1.NumMonomials);
    localSpace.QmatrixInvKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                                  reference_element_data.MonomialsKp1.NumMonomials);




    localSpace.TkNabla.resize(localSpace.NkNabla ,
                              localSpace.Dimension * localSpace.Nk);

    const MatrixXd dx = (1.0 / polygonDiameter) * monomials.D_x(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla,
                                                                                                                      localSpace.Nk);
    const MatrixXd dy = (1.0 / polygonDiameter) * monomials.D_y(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla,
                                                                                                                      localSpace.Nk);

    localSpace.TkNabla << dx , dy;

    MatrixXd V;
    VectorXd S;
    LAPACK_utilities::svd(localSpace.TkNabla, V, S);

    MatrixXd VanderInternal2k(localSpace.Dimension * localSpace.Nk,
                              localSpace.Dimension * internalQuadratureWeights.size());
    VanderInternal2k << localSpace.VanderInternal.transpose(), MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()),
        MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows()), localSpace.VanderInternal.transpose();

    const MatrixXd GkNablaVanderInternal = localSpace.TkNabla * VanderInternal2k;

    localSpace.TkBigOPlus = V.transpose().rightCols(localSpace.Dimension * localSpace.Nk - localSpace.NkNabla).transpose(); // posso usare la base rat?

    const MatrixXd GkBigOPlusVanderInternal = localSpace.TkBigOPlus * VanderInternal2k;

    localSpace.GkVanderInternal.resize(localSpace.Dimension * localSpace.Nk,
                                       localSpace.Dimension * internalQuadratureWeights.size());

    localSpace.GkVanderInternal << GkNablaVanderInternal,
        GkBigOPlusVanderInternal;

    VectorXd internalWeights2k(localSpace.Dimension * internalQuadratureWeights.size());
    internalWeights2k << sqrtInternalQuadratureWeights, sqrtInternalQuadratureWeights;

    const MatrixXd tempG = internalWeights2k.asDiagonal() * localSpace.GkVanderInternal.transpose();
    localSpace.Gmatrix = tempG.transpose() * tempG;
}
//****************************************************************************
void VEM_MCC_2D_VelocityLocalSpace::ComputeDivergenceCoefficients(const double& polytopeMeasure,
                                                                  const Eigen::MatrixXd& W2,
                                                                  VEM_MCC_VelocityLocalSpace_Data& localSpace) const
{
    MatrixXd W1 = MatrixXd::Zero(localSpace.Nk,
                                 localSpace.NumBasisFunctions);

    W1.block(1,
             localSpace.NumBoundaryBasisFunctions,
             localSpace.NumNablaInternalBasisFunctions,
             localSpace.NumNablaInternalBasisFunctions)
        = - polytopeMeasure * Eigen::MatrixXd::Identity(localSpace.NumNablaInternalBasisFunctions,
                                                       localSpace.NumNablaInternalBasisFunctions);


    localSpace.Wmatrix = W1 + W2;
    localSpace.VMatrix = localSpace.Hmatrix.llt().solve(localSpace.Wmatrix);
}
//****************************************************************************
void VEM_MCC_2D_VelocityLocalSpace::ComputeL2Projectors(const double& polygonMeasure,
                                                        const Eigen::MatrixXd& B2Nabla,
                                                        VEM_MCC_VelocityLocalSpace_Data& localSpace) const
{    
    const MatrixXd VanderInternalkp1 = monomialskp1.Vander(internalQuadraturePoints,
                                                           polytopeCentroid,
                                                           polytopeDiameter);

    const MatrixXd HHashtagMatrix = VanderInternalkp1.rightCols(localSpace.NkNabla).transpose() *
                                    internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    const MatrixXd B1Nabla = - HHashtagMatrix * localSpace.VMatrix;

    const MatrixXd VanderBoundarykp1 = monomialskp1.Vander(boundaryQuadraturePoints,
                                                           polytopeCentroid,
                                                           polytopeDiameter);

    const MatrixXd BNabla = B1Nabla + B2Nabla;

    MatrixXd BBigOPlus = MatrixXd::Zero(localSpace.NumBigOPlusInternalBasisFunctions,localSpace.NumBasisFunctions);
    BBigOPlus.rightCols(localSpace.NumBigOPlusInternalBasisFunctions)
        = polygonMeasure * Eigen::MatrixXd::Identity(localSpace.NumBigOPlusInternalBasisFunctions,localSpace.NumBigOPlusInternalBasisFunctions);

    localSpace.Bmatrix.resize(localSpace.Dimension * localSpace.Nk, localSpace.NumBasisFunctions);
    localSpace.Bmatrix << BNabla,
        BBigOPlus;

    localSpace.Pi0k = localSpace.Gmatrix.llt().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_MCC_2D_VelocityLocalSpace::ComputeValuesOnBoundary(const Eigen::MatrixXd& polytopeVertices,
                                                            const Eigen::MatrixXd& edgeNormals,
                                                            const std::vector<bool>& edgeDirections,
                                                            const Eigen::VectorXd& boundaryQuadratureWeights,
                                                            Eigen::VectorXd& edgeDirectionsVector,
                                                            MatrixXd& W2,
                                                            MatrixXd& B2Nabla,
                                                            VEM_MCC_VelocityLocalSpace_Data& localSpace) const
{

    const unsigned int numVertices = polytopeVertices.cols();
    const unsigned int numEdges = numVertices;

    vector<Eigen::VectorXd> edgeNormalsVector(2, Eigen::VectorXd::Zero(localSpace.NumEdgeBasisFunctions));
    edgeDirectionsVector.resize(localSpace.NumEdgeBasisFunctions);

    // offset used below to set edge-internal quadrature points and weights.
    unsigned int edgeInternalPointsOffset = 0;
    for(unsigned int i = 0; i < numEdges; ++i)
    {
        const Eigen::VectorXd& outNormalTimesAbsMapDeterminant = edgeNormals.col(i);
        const double direction = edgeDirections[i] ? 1.0 : -1.0;

        // map edge internal quadrature points
        const unsigned int numEdgeInternalQuadraturePoints = referenceElementVelocity.Order() + 1;

        for(unsigned int d = 0; d < 2; ++d)
        {
            edgeNormalsVector[d].segment(edgeInternalPointsOffset,
                                         numEdgeInternalQuadraturePoints) =
                Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints,
                                          outNormalTimesAbsMapDeterminant[d] * direction);

        }
        edgeDirectionsVector.segment(edgeInternalPointsOffset,
                                     numEdgeInternalQuadraturePoints) =
            Eigen::VectorXd::Constant(numEdgeInternalQuadraturePoints,
                                      direction);

        edgeInternalPointsOffset += numEdgeInternalQuadraturePoints;
    }

    B2Nabla = MatrixXd::Zero(localSpace.NkNabla,localSpace.NumBasisFunctions);
    B2Nabla.block(0,0,
                  localSpace.NkNabla,localSpace.NumEdgeBasisFunctions) =  VanderBoundarykp1.rightCols(localSpace.NkNabla).transpose() *
          boundaryQuadratureWeights.cwiseProduct(edgeDirectionsVector).asDiagonal();

    MatrixXd concatenateEdgeNormalMatrix(localSpace.Dimension * localSpace.NumEdgeBasisFunctions, localSpace.NumEdgeBasisFunctions);
    const MatrixXd temp1 = edgeNormalsVector[0].asDiagonal();
    const MatrixXd temp2 = edgeNormalsVector[1].asDiagonal();
    concatenateEdgeNormalMatrix << temp1,
        temp2;

    MatrixXd VanderBoundary2k(localSpace.Dimension * NumberPolynomialBasisFunctions(),
                              localSpace.Dimension * boundaryQuadratureWeights.size());

    VanderBoundary2k << localSpace.VanderBoundary.transpose(), MatrixXd::Zero(localSpace.VanderBoundary.cols(), localSpace.VanderBoundary.rows()),
        MatrixXd::Zero(localSpace.VanderBoundary.cols(), localSpace.VanderBoundary.rows()), localSpace.VanderBoundary.transpose();

    const MatrixXd GkNablaVanderBoundary = localSpace.TkNabla * VanderBoundary2k;

    const MatrixXd GkBigOPlusVanderBoundary = localSpace.TkBigOPlus * VanderBoundary2k;

    MatrixXd GkVanderBoundary(localSpace.Dimension * NumberPolynomialBasisFunctions(),
                              localSpace.Dimension * boundaryQuadratureWeights.size());

    GkVanderBoundary << GkNablaVanderBoundary,
        GkBigOPlusVanderBoundary;

    localSpace.GkVanderBoundaryTimesNormal = GkVanderBoundary * concatenateEdgeNormalMatrix;
}
//****************************************************************************
}
}
}
