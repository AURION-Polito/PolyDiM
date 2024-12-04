#include "VEM_PCC_2D_Ortho_LocalSpace.hpp"


#include "LAPACK_utilities.hpp"

using namespace Eigen;

namespace Polydim
{
  namespace VEM
  {
    namespace PCC
    {
      //****************************************************************************
      VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Ortho_LocalSpace::CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                               const VEM_PCC_2D_Polygon_Geometry& polygon) const
      {
        VEM_PCC_2D_LocalSpace_Data localSpace;

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

        localSpace.PiNabla = ComputePiNabla(reference_element_data,
                                            polygon.Measure,
                                            polygon.Diameter,
                                            localSpace.InternalQuadrature.Weights,
                                            localSpace.BoundaryQuadrature.Quadrature.Weights,
                                            localSpace.BoundaryQuadrature.WeightsTimesNormal,
                                            localSpace);
        localSpace.StabMatrix = ComputeStabilizationMatrix(polygon.Measure,
                                                           polygon.Diameter,
                                                           localSpace);

        ComputeL2Projectors(polygon.Measure,
                            localSpace);

        ComputeL2ProjectorsOfDerivatives(polygon.Measure,
                                         polygon.Diameter,
                                         localSpace.BoundaryQuadrature.WeightsTimesNormal,
                                         localSpace);

        localSpace.StabMatrixPi0k = ComputeStabilizationMatrixPi0k(polygon.Measure,
                                                                   polygon.Diameter,
                                                                   localSpace);

        return localSpace;
      }
      //****************************************************************************
      VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_Ortho_LocalSpace::Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                                 const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                                                 const Eigen::MatrixXd& internalQuadraturePoints,
                                                                                 const Eigen::VectorXd& internalQuadratureWeights,
                                                                                 const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                                                 const Eigen::VectorXd& boundaryQuadratureWeights,
                                                                                 const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal) const
      {
        VEM_PCC_2D_LocalSpace_Data localSpace;

        InitializeProjectorsComputation(reference_element_data,
                                        polygon.Vertices,
                                        polygon.Centroid,
                                        polygon.Diameter,
                                        internalQuadraturePoints,
                                        internalQuadratureWeights,
                                        boundaryQuadraturePoints,
                                        localSpace);

        localSpace.PiNabla = ComputePiNabla(reference_element_data,
                                            polygon.Measure,
                                            polygon.Diameter,
                                            internalQuadratureWeights,
                                            boundaryQuadratureWeights,
                                            boundaryQuadratureWeightsTimesNormal,
                                            localSpace);

        ComputeL2Projectors(polygon.Measure,
                            localSpace);

        return localSpace;
      }
      //****************************************************************************
      void VEM_PCC_2D_Ortho_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                        const Eigen::MatrixXd& polygonVertices,
                                                                        const Eigen::Vector3d& polygonCentroid,
                                                                        const double& polygonDiameter,
                                                                        const Eigen::MatrixXd& internalQuadraturePoints,
                                                                        const Eigen::VectorXd& internalQuadratureWeights,
                                                                        const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                                        VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        const unsigned int numVertices = polygonVertices.cols();
        const unsigned int numEdges = numVertices;

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
        localSpace.VanderInternal = monomials.Vander(internalQuadraturePoints,
                                                     polygonCentroid,
                                                     polygonDiameter);
        localSpace.VanderInternalDerivatives =  monomials.VanderDerivatives(localSpace.VanderInternal,
                                                                            polygonDiameter);

        localSpace.VanderBoundary = monomials.Vander(boundaryQuadraturePoints,
                                                     polygonCentroid,
                                                     polygonDiameter);

        localSpace.VanderBoundaryDerivatives = monomials.VanderDerivatives(localSpace.VanderBoundary,
                                                                           polygonDiameter);

        //localSpace.internalQuadratureWeightsSqrt = internalQuadratureWeights.array().sqrt();

        // Compute mass matrix of polynomials.
        ChangeOfBasis(internalQuadratureWeights,
                      localSpace);

      }
      //****************************************************************************
      void VEM_PCC_2D_Ortho_LocalSpace::ChangeOfBasis(const Eigen::VectorXd& internalQuadratureWeights,
                                                      VEM_PCC_2D_LocalSpace_Data& localSpace) const
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
      Eigen::MatrixXd VEM_PCC_2D_Ortho_LocalSpace::ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                  const double& polygonMeasure,
                                                                  const double& polygonDiameter,
                                                                  const Eigen::VectorXd& internalQuadratureWeights,
                                                                  const Eigen::VectorXd& boundaryQuadratureWeights,
                                                                  const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                                  VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {

        // G_{ij} = \int_E \nabla m_i \nabla m_j
        MatrixXd Gmatrix = localSpace.VanderInternalDerivatives[0].transpose()*internalQuadratureWeights.asDiagonal() *
                           localSpace.VanderInternalDerivatives[0] + localSpace.VanderInternalDerivatives[1].transpose() *
                           internalQuadratureWeights.asDiagonal()*localSpace.VanderInternalDerivatives[1];
        // B_{ij} = \int_E \nabla m_i \nabla \phi_j
        MatrixXd Bmatrix(localSpace.NumProjectorBasisFunctions,
                         localSpace.NumBasisFunctions);
        // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
        Bmatrix.leftCols(localSpace.NumVertexBasisFunctions +
                         localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundaryDerivatives[0].transpose()*boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
            localSpace.VanderBoundaryDerivatives[1].transpose()*boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

        if(reference_element_data.Order == 1)
        {
          // B_{0j} = \int_{\partial E} \phi_j
          Bmatrix.row(0) = boundaryQuadratureWeights;
          // G_{0j} = \int_{\partial E} m_j
          Gmatrix.row(0) = localSpace.VanderBoundary.transpose()
                           * boundaryQuadratureWeights;
        }
        else
        {
          // G_{0j} = \int_{E} m_j
          Gmatrix.row(0) = localSpace.VanderInternal.transpose()*internalQuadratureWeights;
          // Second block of B: - \int_E \Delta m_i \phi_j
          Bmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
              (- polygonMeasure / (polygonDiameter * polygonDiameter)) *
              (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions)) *
              localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2,
                                                  localSpace.Nkm2);
          // B_{0j} = \int_{E} \phi_j (only the first internal basis
          // function has a non-zero integral)
          //Bmatrix(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) = measure * Qmatrix.inverse().topRows(1).sum();
          MatrixXd a = localSpace.QmatrixInv.row(0);
          Bmatrix.row(0) << MatrixXd::Zero(1,
                                           localSpace.NumVertexBasisFunctions +
                                           localSpace.NumEdgeBasisFunctions),
              polygonMeasure * a.leftCols(localSpace.Nkm2);
        }

        MatrixXd pGmatrix;
        MatrixXd pBmatrix;

        //pGmatrix =Qmatrix*pGmatrix*Qmatrix.transpose();
        MatrixXd temp1 = internalQuadratureWeights.array().sqrt().matrix().asDiagonal() *
                         localSpace.VanderInternalDerivatives[0] *
                         localSpace.Qmatrix.transpose();
        MatrixXd temp2 = internalQuadratureWeights.array().sqrt().matrix().asDiagonal() *
                         localSpace.VanderInternalDerivatives[1] *
                         localSpace.Qmatrix.transpose();
        pGmatrix = temp1.transpose() * temp1 +
                   temp2.transpose() * temp2;

        // Effetuo la QR-rank-revealing fattorizzazione della matrice pTildeGmatrix : pTildeGmatrix * P= Q * R
        ColPivHouseholderQR<MatrixXd> dec(pGmatrix.transpose());
        if (dec.info() != Success) abort();
        MatrixXd P = dec.colsPermutation();
        MatrixXd prod =P.transpose()*pGmatrix;
        prod.bottomLeftCorner(1,
                              localSpace.NumProjectorBasisFunctions) = Gmatrix.topLeftCorner(1,
                                                                                             localSpace.NumProjectorBasisFunctions)
                                                                       * localSpace.Qmatrix.transpose();
        pGmatrix=P*prod;

        pBmatrix = MatrixXd::Zero(localSpace.NumProjectorBasisFunctions,
                                  localSpace.NumBasisFunctions);
        pBmatrix.bottomRightCorner(localSpace.NumProjectorBasisFunctions - 1,
                                   localSpace.NumBasisFunctions) =
            Bmatrix.bottomRightCorner(localSpace.NumProjectorBasisFunctions - 1,
                                      localSpace.NumBasisFunctions);
        pBmatrix = localSpace.Qmatrix * pBmatrix;

        MatrixXd prod1 = P.transpose() * pBmatrix;
        prod1.bottomLeftCorner(1,
                               localSpace.NumBasisFunctions) = Bmatrix.topLeftCorner(1,
                                                                                     localSpace.NumBasisFunctions);
        pBmatrix = P * prod1;

        localSpace.Bmatrix = pBmatrix;
        localSpace.Gmatrix = pGmatrix;

        return pGmatrix.partialPivLu().solve(pBmatrix);
      }
      //****************************************************************************
      MatrixXd VEM_PCC_2D_Ortho_LocalSpace::ComputePolynomialsDofs(const double& polytopeMeasure,
                                                                   const VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        MatrixXd polynomialBasisDofs;
        polynomialBasisDofs.resize(localSpace.NumBasisFunctions,
                                   localSpace.NumProjectorBasisFunctions);
        // boundary degrees of freedom of monomials (values at points on
        // the boundary)
        polynomialBasisDofs.topRows(localSpace.NumVertexBasisFunctions +
                                    localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary * localSpace.Qmatrix.transpose();

        if (referenceElement.Order() > 1)
        {
          // internal degrees of freedom of monomials (scaled moments)
          polynomialBasisDofs.bottomRows(localSpace.NumInternalBasisFunctions) =
              localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polytopeMeasure;
        }

        return polynomialBasisDofs;
      }
      //****************************************************************************
      void VEM_PCC_2D_Ortho_LocalSpace::ComputeL2Projectors(const double& polygonMeasure,
                                                            VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        Eigen::MatrixXd Cmatrix(localSpace.NumProjectorBasisFunctions,
                                localSpace.NumBasisFunctions);
        // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
        Cmatrix.bottomRows(localSpace.NumProjectorBasisFunctions -
                           localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.bottomRows(localSpace.NumProjectorBasisFunctions -
                                          localSpace.NumInternalBasisFunctions) *
            localSpace.PiNabla;

        if(referenceElement.Order() > 1)
        {
          Cmatrix.topLeftCorner(localSpace.NumInternalBasisFunctions,
                                localSpace.NumBasisFunctions -
                                localSpace.NumInternalBasisFunctions).setZero();
          //\int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
          Cmatrix.topRightCorner(localSpace.NumInternalBasisFunctions,
                                 localSpace.NumInternalBasisFunctions) =
              polygonMeasure * Eigen::MatrixXd::Identity(localSpace.NumInternalBasisFunctions,
                                                         localSpace.NumInternalBasisFunctions);
        }

        localSpace.Pi0km1 = localSpace.H_km1_LLT.solve(Cmatrix.topRows(localSpace.Nkm1));
        localSpace.Pi0k = localSpace.Hmatrix.llt().solve(Cmatrix);

        localSpace.Cmatrix = Cmatrix;
      }
      //****************************************************************************
      void VEM_PCC_2D_Ortho_LocalSpace::ComputeL2ProjectorsOfDerivatives(const double& polygonMeasure,
                                                                         const double& polygonDiameter,
                                                                         const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                                         VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        MatrixXd EXmatrix(localSpace.Nkm1,
                          localSpace.NumBasisFunctions);
        EXmatrix.leftCols(localSpace.NumVertexBasisFunctions +
                          localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
            boundaryQuadratureWeightsTimesNormal[0].asDiagonal();
        MatrixXd EYmatrix(localSpace.Nkm1,
                          localSpace.NumBasisFunctions);
        EYmatrix.leftCols(localSpace.NumVertexBasisFunctions +
                          localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
            boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

        if(referenceElement.Order() > 1)
        {
          EXmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
              -(polygonMeasure / polygonDiameter) * monomials.D_x().topLeftCorner(localSpace.Nkm1,
                                                                                  localSpace.NumInternalBasisFunctions) *
              localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2,
                                                  localSpace.Nkm2);
          EYmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
              -(polygonMeasure / polygonDiameter) * monomials.D_y().topLeftCorner(localSpace.Nkm1,
                                                                                  localSpace.NumInternalBasisFunctions) *
              localSpace.QmatrixInv.topLeftCorner(localSpace.Nkm2,
                                                  localSpace.Nkm2);
        }

        localSpace.Pi0km1Der.resize(2);
        localSpace.Pi0km1Der[0] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                   localSpace.Nkm1) * EXmatrix;
        localSpace.Pi0km1Der[1] = localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                   localSpace.Nkm1) * EYmatrix;

        localSpace.Ematrix.resize(2, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));
        localSpace.Ematrix[0] = localSpace.Pi0km1Der[0];
        localSpace.Ematrix[1] = localSpace.Pi0km1Der[1];

        localSpace.Pi0km1Der[0] = localSpace.H_km1_LLT.solve(localSpace.Pi0km1Der[0]);
        localSpace.Pi0km1Der[1] = localSpace.H_km1_LLT.solve(localSpace.Pi0km1Der[1]);
      }

    }
  }
}
