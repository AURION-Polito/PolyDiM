#include "VEM_PCC_3D_Ortho_LocalSpace.hpp"


#include "LAPACK_utilities.hpp"

using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace PCC
{
//****************************************************************************
VEM_PCC_3D_LocalSpace_Data VEM_PCC_3D_Ortho_LocalSpace::CreateLocalSpace(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                         const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    VEM_PCC_3D_LocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature2D;
    const unsigned int numFaces = polyhedron.Faces.size();

    std::vector<Gedim::Quadrature::QuadratureData> FacesInternalQuadrature(numFaces);
    std::vector<Quadrature::VEM_Quadrature_2D::Edges_QuadratureData> FacesBoundaryQuadrature(numFaces);

    for(unsigned int f = 0; f < numFaces; f++)
    {

        localSpace.FacesInternalQuadrature[f]
            = quadrature2D.PolygonInternalQuadrature(reference_element_data.Quadrature,
                                                     polygon.TriangulationVertices);

        localSpace.FacesBoundaryQuadrature[f]
            = quadrature2D.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature,
                                                         polygon.Vertices,
                                                         polygon.EdgesLength,
                                                         polygon.EdgesDirection,
                                                         polygon.EdgesTangent,
                                                         polygon.EdgesNormal);
    }

    Quadrature::VEM_Quadrature_3D quadrature3D;
    localSpace.InternalQuadrature = quadrature3D.PolyhedronInternalQuadrature(reference_element_data.Quadrature,
                                                                              polyhedron.TriangulationVertices);

    localSpace.BoundaryQuadrature = quadrature3D.polyhedronEdgesLobattoQuadrature(reference_element_data.Quadrature,
                                                                                  polyhedron.Vertices,
                                                                                  polyhedron.EdgesLength,
                                                                                  polyhedron.EdgesDirection,
                                                                                  polyhedron.EdgesTangent,
                                                                                  polyhedron.EdgesNormal);


    InitializeProjectorsComputation(reference_element_data,
                                    polyhedron.Vertices,
                                    polyhedron.Edges,
                                    polyhedron.Faces,
                                    polyhedron.Centroid,
                                    polyhedron.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    ComputePiNabla(reference_element_data,
                   polyhedron.Measure,
                   polyhedron.Diameter,
                   localSpace.InternalQuadrature.Weights,
                   localSpace.BoundaryQuadrature.Quadrature.Weights,
                   localSpace.BoundaryQuadrature.WeightsTimesNormal,
                   localSpace);

    ComputeStabilizationMatrix(polyhedron.Diameter,
                               localSpace);

    ComputeL2Projectors(polyhedron.Measure,
                        localSpace);

    ComputeL2ProjectorsOfDerivatives(reference_element_data,
                                     polyhedron.Measure,
                                     polyhedron.Diameter,
                                     localSpace.BoundaryQuadrature.WeightsTimesNormal,
                                     localSpace);

    ComputePolynomialsDofs(polyhedron.Measure,
                           localSpace);

    ComputeStabilizationMatrix(polyhedron.Diameter,
                               localSpace);

    ComputeStabilizationMatrixPi0k(polyhedron.Measure,
                                   localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_PCC_3D_Ortho_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                  const Eigen::MatrixXd& polyhedronVertices,
                                                                  const Eigen::Vector3d& polyhedronCentroid,
                                                                  const double& polyhedronDiameter,
                                                                  const Eigen::MatrixXd& internalQuadraturePoints,
                                                                  const Eigen::VectorXd& internalQuadratureWeights,
                                                                  const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                                  VEM_PCC_3D_LocalSpace_Data& localSpace) const
{
    const unsigned int numVertices = polyhedronVertices.cols();
    const unsigned int numEdges = numVertices;

    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumVertexBasisFunctions = numVertices;
    localSpace.NumEdgeBasisFunctions = (reference_element_data.Order - 1) * numEdges;
    localSpace.NumInternalBasisFunctions = reference_element_data.Order *
                                           (reference_element_data.Order - 1) / 2;

    localSpace.NumBasisFunctions = localSpace.NumVertexBasisFunctions +
                                   localSpace.NumEdgeBasisFunctions +
                                   localSpace.NumInternalBasisFunctions;
    localSpace.NumProjectorBasisFunctions = reference_element_data.Monomials.NumMonomials;

    localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - reference_element_data.Order - 1;
    localSpace.Nkm2 = (reference_element_data.Order - 1) *
                      reference_element_data.Order / 2;

    localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions +
                                           localSpace.NumEdgeBasisFunctions;

    // Compute Vandermonde matrices.
    localSpace.VanderInternal = monomials.Vander(reference_element_data.Monomials,
                                                 internalQuadraturePoints,
                                                 polyhedronCentroid,
                                                 polyhedronDiameter);
    localSpace.VanderInternalDerivatives =  monomials.VanderDerivatives(reference_element_data.Monomials,
                                                                       localSpace.VanderInternal,
                                                                       polyhedronDiameter);

    localSpace.VanderBoundary = monomials.Vander(reference_element_data.Monomials,
                                                 boundaryQuadraturePoints,
                                                 polyhedronCentroid,
                                                 polyhedronDiameter);

    localSpace.VanderBoundaryDerivatives = monomials.VanderDerivatives(reference_element_data.Monomials,
                                                                       localSpace.VanderBoundary,
                                                                       polyhedronDiameter);

    //localSpace.internalQuadratureWeightsSqrt = internalQuadratureWeights.array().sqrt();

    // Compute mass matrix of polynomials.
    ChangeOfBasis(internalQuadratureWeights,
                  localSpace);

}
//****************************************************************************
void VEM_PCC_3D_Ortho_LocalSpace::ChangeOfBasis(const Eigen::VectorXd& internalQuadratureWeights,
                                                VEM_PCC_3D_LocalSpace_Data& localSpace) const
{
    MatrixXd Q1;
    MatrixXd R1;
    LAPACK_utilities::MGS(localSpace.VanderInternal,
                          Q1,
                          R1);

    // L2(E)-re-orthogonalization process
    MatrixXd Q2;
    MatrixXd R2;
    LAPACK_utilities::MGS(internalQuadratureWeights.array().sqrt().matrix().asDiagonal() * Q1,
                          Q2,
                          R2);

    localSpace.Hmatrix= Q2.transpose() * Q2;
    localSpace.H_km1_LLT =  localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1,
                                                            localSpace.Nkm1).llt();

    localSpace.QmatrixInv = (R2 * R1).transpose();
    LAPACK_utilities::inverseTri(localSpace.QmatrixInv,
                                 localSpace.Qmatrix,
                                 'L', 'N');

}
//****************************************************************************
void VEM_PCC_3D_Ortho_LocalSpace::ComputePiNabla(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                 const double& polyhedronMeasure,
                                                 const double& polyhedronDiameter,
                                                 const Eigen::VectorXd& internalQuadratureWeights,
                                                 const Eigen::VectorXd& boundaryQuadratureWeights,
                                                 const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                 VEM_PCC_3D_LocalSpace_Data& localSpace) const
{

    // G_{ij} = \int_E \nabla m_i \nabla m_j
    localSpace.Gmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions,localSpace.NumProjectorBasisFunctions);

    const VectorXd internalQuadratureWeightsSqrt = internalQuadratureWeights.array().sqrt();
    for(unsigned int d = 0; d < localSpace.Dimension; d++)
    {
        const MatrixXd temp = internalQuadratureWeightsSqrt.asDiagonal()
                              * localSpace.VanderInternalDerivatives[d]
                              * localSpace.Qmatrix.transpose();

        localSpace.Gmatrix += temp.transpose() * temp;
    }

    // B_{ij} = \int_E \nabla m_i \nabla \phi_j
    localSpace.Bmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions,
                                        localSpace.NumBasisFunctions);
    // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
    localSpace.Bmatrix.leftCols(localSpace.NumBoundaryBasisFunctions) =
        (localSpace.VanderBoundaryDerivatives[0].transpose()*boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
         localSpace.VanderBoundaryDerivatives[1].transpose()*boundaryQuadratureWeightsTimesNormal[1].asDiagonal() +
         localSpace.VanderBoundaryDerivatives[2].transpose()*boundaryQuadratureWeightsTimesNormal[2].asDiagonal()) *
        localSpace.VanderFaceProjections;

    if(localSpace.Order == 1)
    {
        // B_{0j} = \int_{\partial E} \phi_j
        localSpace.Bmatrix.row(0) = boundaryQuadratureWeights.transpose() *
                                    localSpace.VanderFaceProjections;;
        // G_{0j} = \int_{\partial E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderBoundary.transpose() * boundaryQuadratureWeights;
    }
    else
    {
        // G_{0j} = \int_{E} m_j
        localSpace.Gmatrix.row(0) = localSpace.VanderInternal.transpose()*internalQuadratureWeights;
        // Second block of B: - \int_E \Delta m_i \phi_j
        localSpace.Bmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
            (- polyhedronMeasure / (polyhedronDiameter * polyhedronDiameter)) *
            (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions)) *
            localSpace.QmatrixInv.topLeftCorner(localSpace.NumInternalBasisFunctions,
                                                localSpace.NumInternalBasisFunctions);
        // B_{0j} = \int_{E} \phi_j (only the first internal basis
        // function has a non-zero integral)
        //Bmatrix(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) = measure * Qmatrix.inverse().topRows(1).sum();
        localSpace.Bmatrix.row(0) << MatrixXd::Zero(1,
                                                    localSpace.NumBoundaryBasisFunctions),
            polyhedronMeasure * localSpace.QmatrixInv.row(0).leftCols(localSpace.NumInternalBasisFunctions);
    }

    localSpace.Bmatrix = localSpace.Qmatrix * localSpace.Bmatrix;

    localSpace.PiNabla =  localSpace.Gmatrix.partialPivLu().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_PCC_3D_Ortho_LocalSpace::ComputePolynomialsDofs(const double& polyhedronMeasure,
                                                         VEM_PCC_3D_LocalSpace_Data& localSpace) const
{
    localSpace.Dmatrix.setZero(localSpace.NumBasisFunctions,
                               localSpace.NumProjectorBasisFunctions);

    localSpace.Dmatrix.topRows(localSpace.VanderEdgeDofs.rows()) = localSpace.VanderEdgeDofs;

    if(localSpace.Order > 1)
    {
        // internal scaled moments on each face.
        localSpace.Dmatrix.block(localSpace.VanderEdgeDofs.rows(), 0,
                                 localSpace.NumFaceBasisFunctions,
                                 localSpace.NumProjectorBasisFunctions) = localSpace.ScaledHmatrixOnBoundary;


        // internal scaled moments
        localSpace.Dmatrix.bottomRows(localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polyhedronMeasure;
    }
}
//****************************************************************************
void VEM_PCC_3D_Ortho_LocalSpace::ComputeL2ProjectorsOfDerivatives(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                   const double& polyhedronMeasure,
                                                                   const double& polyhedronDiameter,
                                                                   const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                                   VEM_PCC_3D_LocalSpace_Data& localSpace) const
{
    localSpace.Pi0km1Der.resize(localSpace.Dimension);

    localSpace.Ematrix.resize(3, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));

    for ( unsigned int d = 0; d < localSpace.Dimension; d++ )
    {
        localSpace.Ematrix[d].leftCols(localSpace.NumBasisFunctions - localSpace.NumInternalBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose()
            * boundaryQuadratureWeightsTimesNormal[d].asDiagonal()
            * localSpace.VanderFaceProjections;

        if(localSpace.Order > 1)
        {
            localSpace.Ematrix[d].rightCols(localSpace.NumInternalBasisFunctions) =
                -(polyhedronMeasure/polyhedronDiameter)
                * monomials.DerivativeMatrix(reference_element_data.Monomials, d).topLeftCorner(localSpace.Nkm1, localSpace.NumInternalBasisFunctions)
                * localSpace.QmatrixInv.topLeftCorner(localSpace.NumInternalBasisFunctions, localSpace.NumInternalBasisFunctions);
        }

        localSpace.Ematrix[d] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1)
                                * localSpace.Ematrix[d];

        localSpace.Pi0km1Der[d] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[d]);

    }
}
//****************************************************************************

}
}
}
