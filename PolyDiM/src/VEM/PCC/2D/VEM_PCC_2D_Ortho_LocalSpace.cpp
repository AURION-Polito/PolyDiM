#include "VEM_PCC_2D_Ortho_LocalSpace.hpp"

#include "LAPACK_utilities.hpp"
#include "iostream"

using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Ortho_LocalSpace::CreateLocalSpace(
    const VEM_PCC_2D_ReferenceElement_Data &reference_element_data, const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(
        reference_element_data.Quadrature.ReferenceTriangleQuadrature, polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
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
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Ortho_LocalSpace::Compute3DUtilities(
    const VEM_PCC_2D_ReferenceElement_Data &reference_element_data, const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(
        reference_element_data.Quadrature.ReferenceTriangleQuadrature, polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceSegmentInternalPoints,
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

    return localSpace;
}
//****************************************************************************
void VEM_PCC_2D_Ortho_LocalSpace::InitializeProjectorsComputation(
    const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
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
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;

    localSpace.NumBasisFunctions =
        localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions + localSpace.NumInternalBasisFunctions;
    localSpace.NumProjectorBasisFunctions = reference_element_data.Monomials.NumMonomials;

    localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - reference_element_data.Order - 1;
    localSpace.Nkm2 = (reference_element_data.Order - 1) * reference_element_data.Order / 2;

    localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions;

    // Compute Vandermonde matrices.
    localSpace.Diameter = polygonDiameter;
    localSpace.Centroid = polygonCentroid;
    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polygonCentroid, polygonDiameter);
    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polygonDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polygonDiameter);

    // Compute mass matrix of polynomials.
    ChangeOfBasis(internalQuadratureWeights, localSpace);
}
//****************************************************************************
void VEM_PCC_2D_Ortho_LocalSpace::ChangeOfBasis(const Eigen::VectorXd &internalQuadratureWeights,
                                                VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    MatrixXd Q1;
    MatrixXd R1;
    LAPACK_utilities::MGS(localSpace.VanderInternal, Q1, R1);

    // L2(E)-re-orthogonalization process
    MatrixXd Q2;
    MatrixXd R2;
    LAPACK_utilities::MGS(internalQuadratureWeights.array().sqrt().matrix().asDiagonal() * Q1, Q2, R2);

    localSpace.Hmatrix = Q2.transpose() * Q2;
    localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();

    localSpace.QmatrixInv = (R2 * R1).transpose();
    LAPACK_utilities::inverseTri(localSpace.QmatrixInv, localSpace.Qmatrix, 'L', 'N');
}
//****************************************************************************
void VEM_PCC_2D_Ortho_LocalSpace::ComputePiNabla(
    const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
    const double &polygonMeasure,
    const double &polygonDiameter,
    const Eigen::VectorXd &internalQuadratureWeights,
    const Eigen::VectorXd &boundaryQuadratureWeights,
    const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
    VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    const MatrixXd temp1 = internalQuadratureWeights.array().sqrt().matrix().asDiagonal() *
                           localSpace.VanderInternalDerivatives[0] * localSpace.Qmatrix.transpose();
    const MatrixXd temp2 = internalQuadratureWeights.array().sqrt().matrix().asDiagonal() *
                           localSpace.VanderInternalDerivatives[1] * localSpace.Qmatrix.transpose();

    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = temp1.transpose() * temp1 + temp2.transpose() * temp2;

    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    localSpace.Bmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);
    // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
    localSpace.Bmatrix.leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundaryDerivatives[0].transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
        localSpace.VanderBoundaryDerivatives[1].transpose() * boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (reference_element_data.Order == 1)
    {
        localSpace.Bmatrix = localSpace.Qmatrix * localSpace.Bmatrix;

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
            (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions)) *
            localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2, localSpace.Nkm2);
        // B_{0j} = \int_{E} \phi_j (only the first internal basis
        // function has a non-zero integral)

        localSpace.Bmatrix = localSpace.Qmatrix * localSpace.Bmatrix;
        localSpace.Bmatrix.row(0) << MatrixXd::Zero(
            1, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions),
            polygonMeasure * localSpace.QmatrixInv.row(0).leftCols(localSpace.Nkm2);
    }
    localSpace.Gmatrix.row(0) = localSpace.Gmatrix.row(0) * localSpace.Qmatrix.transpose();

    localSpace.PiNabla = localSpace.Gmatrix.partialPivLu().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_PCC_2D_Ortho_LocalSpace::ComputePolynomialsDofs(const double &polytopeMeasure,
                                                         VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    localSpace.Dmatrix.resize(localSpace.NumBasisFunctions, localSpace.NumProjectorBasisFunctions);
    // boundary degrees of freedom of monomials (values at points on
    // the boundary)
    localSpace.Dmatrix.topRows(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundary * localSpace.Qmatrix.transpose();

    if (localSpace.Order > 1)
    {
        // internal degrees of freedom of monomials (scaled moments)
        localSpace.Dmatrix.bottomRows(localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polytopeMeasure;
    }
}
//****************************************************************************
void VEM_PCC_2D_Ortho_LocalSpace::ComputeL2ProjectorsOfDerivatives(
    const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
    const double &polygonMeasure,
    const double &polygonDiameter,
    const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
    VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    localSpace.Ematrix.resize(2, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));

    localSpace.Ematrix[0].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
        boundaryQuadratureWeightsTimesNormal[0].asDiagonal();

    localSpace.Ematrix[1].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
        localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
        boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

    if (localSpace.Order > 1)
    {
        localSpace.Ematrix[0].rightCols(localSpace.NumInternalBasisFunctions) =
            -(polygonMeasure / polygonDiameter) *
            monomials.D_x(reference_element_data.Monomials)
                .topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions) *
            localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2, localSpace.Nkm2);
        localSpace.Ematrix[1].rightCols(localSpace.NumInternalBasisFunctions) =
            -(polygonMeasure / polygonDiameter) *
            monomials.D_y(reference_element_data.Monomials)
                .topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions) *
            localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2, localSpace.Nkm2);
    }

    localSpace.Ematrix[0] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1) * localSpace.Ematrix[0];
    localSpace.Ematrix[1] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1) * localSpace.Ematrix[1];

    localSpace.Pi0km1Der.resize(2);
    localSpace.Pi0km1Der[0] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[0]);
    localSpace.Pi0km1Der[1] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[1]);
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
