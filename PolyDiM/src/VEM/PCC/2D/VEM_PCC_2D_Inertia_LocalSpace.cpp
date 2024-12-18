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
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometry_utilities(geometryUtilitiesConfig);

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature = quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
                                                                             reference_element_data.Quadrature.ReferenceSegmentInternalWeights,
                                                                             reference_element_data.Quadrature.ReferenceSegmentExtremaWeights,
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

    ComputePiNabla(reference_element_data,
                   polygon.Measure,
                   polygon.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.BoundaryQuadrature.Quadrature.Weights,
                   localSpace.BoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(polygon.Measure, localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data,
                                     polygon.Measure,
                                     polygon.Diameter,
                                     localSpace.BoundaryQuadrature.WeightsTimesNormal,
                                     localSpace);

    ComputePolynomialsDofs(polygon.Measure, localSpace);

    ComputeStabilizationMatrix(polygon.Diameter, localSpace);

    ComputeStabilizationMatrixPi0k(polygon.Measure, localSpace);

    return localSpace;
}
// ***************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::InertiaMapping(const Gedim::GeometryUtilities& geometryUtilities,
                                                   const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                   InertiaData &inertia_data)
{

    inertia_data.FmatrixInv = Matrix3d::Identity();
    inertia_data.Fmatrix = Matrix3d::Identity();
    inertia_data.absDetFmatrix  = 1.0;
    inertia_data.translation = Vector3d::Zero();

    // First Rescaling
    const double invPolygonDiameter = (1.0 / polygon.Diameter);

    inertia_data.FmatrixInv = invPolygonDiameter * inertia_data.FmatrixInv;
    inertia_data.Fmatrix = inertia_data.Fmatrix * polygon.Diameter;
    inertia_data.absDetFmatrix  *= polygon.Diameter *
                     polygon.Diameter;
    inertia_data.translation += Vector3d::Zero();

    const MatrixXd polygonVerticesFirstRescaling = inertia_data.FmatrixInv * (polygon.Vertices.colwise() - inertia_data.translation);
    const Vector3d polygonCentroidsFirstRescaling = inertia_data.FmatrixInv * (polygon.Centroid - inertia_data.translation);

    // Inertia Mapping
    vector<Matrix3d> polygonTriangulationsFirstRescaling;
    unsigned int numTriangles = polygon.TriangulationVertices.size();
    polygonTriangulationsFirstRescaling.resize(numTriangles);
    for(unsigned int n = 0; n < numTriangles; n++)
        polygonTriangulationsFirstRescaling[n] = inertia_data.FmatrixInv * (polygon.TriangulationVertices[n].colwise() - inertia_data.translation);

    const Matrix2d HmatrixFirstRescaling = geometryUtilities.PolygonMass(polygonCentroidsFirstRescaling,
                                                                         polygonTriangulationsFirstRescaling);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(HmatrixFirstRescaling);
    if (eigensolver.info() != Eigen::Success) abort();
    Vector2d sqrtLambdaFirstRescaling =  eigensolver.eigenvalues().array().sqrt();
    Vector2d sqrtLambdaInvFirstRescaling =  eigensolver.eigenvalues().array().rsqrt();
    MatrixXd QmatrixFirstRescaling = eigensolver.eigenvectors();
    MatrixXd Bmatrix = sqrtLambdaInvFirstRescaling.asDiagonal() * QmatrixFirstRescaling.transpose() * sqrtLambdaFirstRescaling.maxCoeff();
    MatrixXd BmatrixInv = QmatrixFirstRescaling * sqrtLambdaFirstRescaling.asDiagonal() / sqrtLambdaFirstRescaling.maxCoeff();

    inertia_data.FmatrixInv.topLeftCorner(2,2) = Bmatrix * inertia_data.FmatrixInv.topLeftCorner(2,2) ;
    inertia_data.Fmatrix.topLeftCorner(2,2) =  inertia_data.Fmatrix.topLeftCorner(2,2) * BmatrixInv;
    const double detB = BmatrixInv.determinant();
    inertia_data.signDetQ = (detB > 1.0e-12) ? 1.0 : -1.0;

#ifdef TEST
    if (abs(detB) <= 1.0e-12)
        throw runtime_error("Singular matrix of inertia");
#endif

    inertia_data.absDetFmatrix  *= abs(detB);
    inertia_data.translation += polygon.Centroid;

    const MatrixXd polygonVerticesInertiaMapping = inertia_data.FmatrixInv * (polygon.Vertices.colwise() - inertia_data.translation);

    // Second rescaling

    const double polygonDiameterInertiaMapping = geometryUtilities.PolygonDiameter(polygonVerticesInertiaMapping);

    const double invpolygonDiameterInertiaMapping = (1.0 / polygonDiameterInertiaMapping);

    inertia_data.FmatrixInv = invpolygonDiameterInertiaMapping * inertia_data.FmatrixInv ;
    inertia_data.Fmatrix = polygonDiameterInertiaMapping * inertia_data.Fmatrix;
    inertia_data.absDetFmatrix  *= polygonDiameterInertiaMapping *
                     polygonDiameterInertiaMapping;
    inertia_data.translation += Vector3d::Zero();

    inertia_data.Vertices = inertia_data.FmatrixInv * (polygon.Vertices.colwise() - inertia_data.translation);

    ComputeGeometryProperties(geometryUtilities,
                              polygon.EdgesDirection,
                              polygon.TriangulationVertices,
                              inertia_data);


}
// ***************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputeGeometryProperties(const Gedim::GeometryUtilities& geometryUtilities,
                                                              const vector<bool>& polygonEdgeDirections,
                                                              const std::vector<Eigen::Matrix3d>& polygonTriangulation,
                                                              InertiaData &data)
{
    // Extract original cell2D geometric information
    unsigned int numVertices = data.Vertices.cols();
    // compute original cell2D triangulation
    if (data.signDetQ < 0)
    {
        data.OrderedVertices.resize(3, numVertices);
        for(unsigned int v = 0; v < numVertices; v++)
            data.OrderedVertices.col(v) = data.Vertices.col(numVertices - v - 1);
    }
    else
        data.OrderedVertices = data.Vertices;

    const unsigned int& cell2DTriangulationSize = polygonTriangulation.size();
    vector<Eigen::Matrix3d> convexCell2DTriangulationPoints(cell2DTriangulationSize);
    for(unsigned int i = 0; i < cell2DTriangulationSize; i++)
    {
        if (data.signDetQ < 0)
            for(unsigned int v = 0; v < 3; v++)
            {
                convexCell2DTriangulationPoints[i].col(v) = data.FmatrixInv * (polygonTriangulation[i].col(2 - v) - data.translation);
            }
        else
            convexCell2DTriangulationPoints[i] = data.FmatrixInv * (polygonTriangulation[i].colwise() - data.translation);
    }

    // compute original cell2D area and centroids
    Eigen::VectorXd convexCell2DTriangulationAreas(cell2DTriangulationSize);
    Eigen::MatrixXd convexCell2DTriangulationCentroids(3, cell2DTriangulationSize);
    for (unsigned int cct = 0; cct < cell2DTriangulationSize; cct++)
    {
        convexCell2DTriangulationAreas[cct] = geometryUtilities.PolygonArea(convexCell2DTriangulationPoints[cct]);
        convexCell2DTriangulationCentroids.col(cct) = geometryUtilities.PolygonCentroid(convexCell2DTriangulationPoints[cct],
                                                                                        convexCell2DTriangulationAreas[cct]);
    }

    data.Measure = convexCell2DTriangulationAreas.sum();
    data.Centroid = geometryUtilities.PolygonCentroid(convexCell2DTriangulationCentroids,
                                                      convexCell2DTriangulationAreas,
                                                      data.Measure);

    // Compute cell2D triangulation from original cell2Ds
    std::vector<Eigen::Matrix3d> mappedPolygonTriangulations(cell2DTriangulationSize);
    unsigned int triangulationCounter = 0;
    for (unsigned int cct = 0; cct < convexCell2DTriangulationPoints.size(); cct++)
        mappedPolygonTriangulations[triangulationCounter++] = convexCell2DTriangulationPoints[cct];

    data.Diameter = geometryUtilities.PolygonDiameter(data.Vertices);

    data.EdgesLength = geometryUtilities.PolygonEdgeLengths(data.Vertices);
    data.EdgesTangent = geometryUtilities.PolygonEdgeTangents(data.Vertices);

    data.EdgesNormal = data.signDetQ *
                       geometryUtilities.PolygonEdgeNormals(data.Vertices);

    data.EdgesDirection = polygonEdgeDirections;

}
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Inertia_LocalSpace::Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                             const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature = quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
                                                                             reference_element_data.Quadrature.ReferenceSegmentInternalWeights,
                                                                             reference_element_data.Quadrature.ReferenceSegmentExtremaWeights,
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

    ComputePiNabla(reference_element_data,
                   polygon.Measure,
                   polygon.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.BoundaryQuadrature.Quadrature.Weights,
                   localSpace.BoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(polygon.Measure,
                        localSpace);

    return localSpace;
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
    localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions
                                           + localSpace.NumEdgeBasisFunctions;
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.NumBasisFunctions = localSpace.NumVertexBasisFunctions
                                   + localSpace.NumEdgeBasisFunctions
                                   + localSpace.NumInternalBasisFunctions;

    localSpace.NumProjectorBasisFunctions = reference_element_data.Monomials.NumMonomials;

    localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - reference_element_data.Order - 1;


    // Compute Vandermonde matrices.
    localSpace.VanderInternal = monomials.Vander(reference_element_data.Monomials,
                                                 internalQuadraturePoints,
                                                 polygonCentroid,
                                                 polygonDiameter);
    localSpace.VanderInternalDerivatives = monomials
                                               .VanderDerivatives(reference_element_data.Monomials,
                                                                  localSpace.VanderInternal,
                                                                  polygonDiameter);

    localSpace.VanderBoundary = monomials.Vander(reference_element_data.Monomials,
                                                 boundaryQuadraturePoints,
                                                 polygonCentroid,
                                                 polygonDiameter);

    localSpace.VanderBoundaryDerivatives = monomials
                                               .VanderDerivatives(reference_element_data.Monomials,
                                                                  localSpace.VanderBoundary,
                                                                  polygonDiameter);

    // Compute mass matrix of monomials.
    localSpace.Hmatrix = localSpace.VanderInternal.transpose()
                         * internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;

    localSpace.Qmatrix = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions,
                                            localSpace.NumProjectorBasisFunctions);
    localSpace.QmatrixInv = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions,
                                               localSpace.NumProjectorBasisFunctions);

    // Compute LLT factorization of order-1 monomials.
    localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();
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
    localSpace.Gmatrix = localSpace.VanderInternalDerivatives[0].transpose()
                             * internalQuadratureWeights.asDiagonal()
                             * localSpace.VanderInternalDerivatives[0]
                         + localSpace.VanderInternalDerivatives[1].transpose()
                               * internalQuadratureWeights.asDiagonal()
                               * localSpace.VanderInternalDerivatives[1];
    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    localSpace.Bmatrix.setZero(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);
    // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
    localSpace.Bmatrix.leftCols(localSpace.NumVertexBasisFunctions
                                + localSpace.NumEdgeBasisFunctions)
        = localSpace.VanderBoundaryDerivatives[0].transpose()
              * boundaryQuadratureWeightsTimesNormal[0].asDiagonal()
          + localSpace.VanderBoundaryDerivatives[1].transpose()
                * boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (localSpace.Order == 1) {
        // B_{0j} = \int_{\partial E} \phi_j
        localSpace.Bmatrix.row(0) = boundaryQuadratureWeights;
        // G_{0j} = \int_{\partial E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderBoundary.transpose()
                                    * boundaryQuadratureWeights;
    } else {
        // G_{0j} = \int_{E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderInternal.transpose()
                                    * internalQuadratureWeights;
        // Second block of B: - \int_E \Delta m_i \phi_j
        localSpace.Bmatrix.rightCols(localSpace.NumInternalBasisFunctions)
            = (-polygonMeasure / (polygonDiameter * polygonDiameter))
              * (reference_element_data.Monomials.Laplacian.leftCols(
                  localSpace.NumInternalBasisFunctions));
        // B_{0j} = \int_{E} \phi_j (only the first internal basis
        // function has a non-zero integral)
        localSpace.Bmatrix(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions)
            = polygonMeasure;
    }

    localSpace.PiNabla = localSpace.Gmatrix.partialPivLu().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_PCC_2D_Inertia_LocalSpace::ComputePolynomialsDofs(const double &polygonMeasure,
                                                           VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix = MatrixXd::Zero(localSpace.NumBasisFunctions,
                                        localSpace.NumProjectorBasisFunctions);
    // boundary degrees of freedom of monomials (values at points on
    // the boundary)
    localSpace.Dmatrix.topRows(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions)
        = localSpace.VanderBoundary;

    if (localSpace.Order > 1) {
        // internal degrees of freedom of monomials (scaled moments)
        localSpace.Dmatrix.bottomRows(localSpace.NumInternalBasisFunctions)
            = localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polygonMeasure;
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
    localSpace.Ematrix[0].leftCols(localSpace.NumVertexBasisFunctions
                                   + localSpace.NumEdgeBasisFunctions)
        = localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose()
          * boundaryQuadratureWeightsTimesNormal[0].asDiagonal();

    localSpace.Ematrix[1].leftCols(localSpace.NumVertexBasisFunctions
                                   + localSpace.NumEdgeBasisFunctions)
        = localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose()
          * boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (localSpace.Order > 1) {
        localSpace.Ematrix[0].rightCols(localSpace.NumInternalBasisFunctions)
            = -(polygonMeasure / polygonDiameter)
              * monomials.D_x(reference_element_data.Monomials)
                    .topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions);
        localSpace.Ematrix[1].rightCols(localSpace.NumInternalBasisFunctions)
            = -(polygonMeasure / polygonDiameter)
              * monomials.D_y(reference_element_data.Monomials)
                    .topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions);
    }

    localSpace.Pi0km1Der.resize(2);
    localSpace.Pi0km1Der[0] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[0]);
    localSpace.Pi0km1Der[1] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[1]);
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
