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

#include "VEM_PCC_2D_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_LocalSpace::CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

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
                                    polygon.Vertices,
                                    polygon.Centroid,
                                    polygon.Measure,
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

    return localSpace;
}
//****************************************************************************
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_LocalSpace::Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

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
                                    polygon.Vertices,
                                    polygon.Centroid,
                                    polygon.Measure,
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
VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_LocalSpace::Compute3DUtilities_DF_PCC(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_PCC_2D_Polygon_Geometry &polygon) const
{
    VEM_PCC_2D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature.ReferenceTriangleQuadrature,
                                                                         polygon.TriangulationVertices);

    localSpace.BoundaryQuadrature =
        quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsInternalWeights,
                                                 reference_element_data.Quadrature.ReferenceEdgeDOFsExtremaWeights,
                                                 polygon.Vertices,
                                                 polygon.EdgesLength,
                                                 polygon.EdgesDirection,
                                                 polygon.EdgesTangent,
                                                 polygon.EdgesNormal);

    Quadrature::VEM_QuadratureData_2D quadrature_data_KL = quadrature.Compute_DF_PCC_3D(reference_element_data.Order);
    localSpace.InternalQuadratureKL =
        quadrature.PolygonInternalQuadrature(quadrature_data_KL.ReferenceTriangleQuadrature, polygon.TriangulationVertices);

    InitializeE2ProjectorsComputation(reference_element_data,
                                      2,
                                      polygon.Vertices,
                                      polygon.Centroid,
                                      polygon.Measure,
                                      polygon.Diameter,
                                      localSpace.InternalQuadrature.Points,
                                      localSpace.InternalQuadrature.Weights,
                                      localSpace.InternalQuadratureKL.Points,
                                      localSpace.InternalQuadratureKL.Weights,
                                      localSpace.BoundaryQuadrature.Quadrature.Points,
                                      localSpace);

    ComputePolynomialsDofs(polygon.Measure, localSpace);

    ComputePiNabla(reference_element_data,
                   polygon.Measure,
                   polygon.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.BoundaryQuadrature.Quadrature.Weights,
                   localSpace.BoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeL2Projectors(polygon.Measure, localSpace);

    ComputeL2ProjectorsKL(localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_PCC_2D_LocalSpace::InitializeE2ProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                              const unsigned int &l,
                                                              const Eigen::MatrixXd &polygonVertices,
                                                              const Eigen::Vector3d &polygonCentroid,
                                                              const double &polygonMeasure,
                                                              const double &polygonDiameter,
                                                              const Eigen::MatrixXd &internalQuadraturePoints,
                                                              const Eigen::VectorXd &internalQuadratureWeights,
                                                              const Eigen::MatrixXd &internalQuadratureKLPoints,
                                                              const Eigen::VectorXd &internalQuadratureKLWeights,
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

    localSpace.Nklm1 = (reference_element_data.Order + l) * (reference_element_data.Order + l + 1) / 2;

    // Compute Vandermonde matrices.
    localSpace.Diameter = polygonDiameter;
    localSpace.Centroid = polygonCentroid;
    localSpace.Measure = polygonMeasure;

    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderInternal, polygonDiameter);

    localSpace.VanderBoundary =
        monomials.Vander(reference_element_data.Monomials, boundaryQuadraturePoints, polygonCentroid, polygonDiameter);

    localSpace.VanderBoundaryDerivatives =
        monomials.VanderDerivatives(reference_element_data.Monomials, localSpace.VanderBoundary, polygonDiameter);

    localSpace.VanderInternalKL =
        monomials.Vander(monomials.Compute(localSpace.Order + l - 1), internalQuadratureKLPoints, polygonCentroid, polygonDiameter);

    // Compute mass matrix of monomials.
    localSpace.Hmatrix = localSpace.VanderInternal.transpose() * internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;

    // Compute LLT factorization of order-1 monomials.
    localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();

    Eigen::MatrixXd H_klm1_matrix;
    if (l == 0)
        localSpace.H_klm1_matrix = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1);
    else
        localSpace.H_klm1_matrix =
            localSpace.VanderInternalKL.transpose() * internalQuadratureKLWeights.asDiagonal() * localSpace.VanderInternalKL;

    // Compute LLT factorization of order-1 monomials.
    localSpace.H_klm1_LLT = localSpace.H_klm1_matrix.llt();
}
//****************************************************************************
void VEM_PCC_2D_LocalSpace::ComputeL2ProjectorsKL(VEM_PCC_2D_LocalSpace_Data &localSpace) const
{
    Eigen::MatrixXd Cmatrix(localSpace.Nklm1, localSpace.NumBasisFunctions);

    Cmatrix.topRows(localSpace.NumProjectorBasisFunctions) = localSpace.Cmatrix;

    Cmatrix.bottomRows(localSpace.Nklm1 - localSpace.NumProjectorBasisFunctions) =
        localSpace.H_klm1_matrix.bottomLeftCorner(localSpace.Nklm1 - localSpace.NumProjectorBasisFunctions,
                                                  localSpace.NumProjectorBasisFunctions) *
        localSpace.PiNabla;

    localSpace.Pi0klm1 = localSpace.H_klm1_LLT.solve(Cmatrix);

    {
        double test_error = (localSpace.H_klm1_matrix.topLeftCorner(localSpace.Nklm1, localSpace.NumProjectorBasisFunctions) -
                             Cmatrix * localSpace.Dmatrix)
                                .norm();
        assert(test_error < 1.0e-12);
    }
}
//****************************************************************************
void VEM_PCC_2D_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                            const Eigen::MatrixXd &polygonVertices,
                                                            const Eigen::Vector3d &polygonCentroid,
                                                            const double &polygonMeasure,
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
    localSpace.Diameter = polygonDiameter;
    localSpace.Centroid = polygonCentroid;
    localSpace.Measure = polygonMeasure;

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
    localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();
}
//****************************************************************************
void VEM_PCC_2D_LocalSpace::ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
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
void VEM_PCC_2D_LocalSpace::ComputePolynomialsDofs(const double &polygonMeasure, VEM_PCC_2D_LocalSpace_Data &localSpace) const
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
void VEM_PCC_2D_LocalSpace::ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
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
    localSpace.Pi0km1Der[0] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[0]);
    localSpace.Pi0km1Der[1] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[1]);
}
//****************************************************************************
} // namespace PCC
} // namespace VEM
} // namespace Polydim
