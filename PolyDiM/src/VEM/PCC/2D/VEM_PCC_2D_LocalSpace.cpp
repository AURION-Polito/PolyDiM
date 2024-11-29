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
      VEM_PCC_2D_LocalSpace::VEM_PCC_2D_LocalSpace()
      {
      }
      //****************************************************************************
      VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_LocalSpace::CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
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

        localSpace.StabMatrix = ComputeStabilizationMatrix(polygon.Measure,
                                                           polygon.Diameter,
                                                           localSpace);

        ComputeL2Projectors(polygon.Measure,
                            localSpace);
        ComputeL2ProjectorsOfDerivatives(reference_element_data,
                                         polygon.Measure,
                                         polygon.Diameter,
                                         boundaryQuadratureWeightsTimesNormal,
                                         localSpace);

        localSpace.StabMatrixPi0k = ComputeStabilizationMatrixPi0k(polygon.Measure,
                                                                   localSpace);

        return localSpace;
      }
      //****************************************************************************
      VEM_PCC_2D_LocalSpace_Data VEM_PCC_2D_LocalSpace::Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
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
      void VEM_PCC_2D_LocalSpace::InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
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
        localSpace.NumInternalBasisFunctions = reference_element_data.Order * (reference_element_data.Order - 1) / 2;

        localSpace.NumBasisFunctions = localSpace.NumVertexBasisFunctions +
                                       localSpace.NumEdgeBasisFunctions +
                                       localSpace.NumInternalBasisFunctions;

        localSpace.NumProjectorBasisFunctions = reference_element_data.Monomials.NumMonomials;

        localSpace.Nkm1 = localSpace.NumProjectorBasisFunctions - reference_element_data.Order - 1;

        localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions +
                                               localSpace.NumEdgeBasisFunctions;

        // Compute Vandermonde matrices.
        localSpace.VanderInternal = monomials.Vander(reference_element_data.Monomials,
                                                     internalQuadraturePoints,
                                                     polygonCentroid,
                                                     polygonDiameter);
        localSpace.VanderInternalDerivatives =  monomials.VanderDerivatives(reference_element_data.Monomials,
                                                                            localSpace.VanderInternal,
                                                                            polygonDiameter);

        localSpace.VanderBoundary = monomials.Vander(reference_element_data.Monomials,
                                                     boundaryQuadraturePoints,
                                                     polygonCentroid,
                                                     polygonDiameter);

        localSpace.VanderBoundaryDerivatives = monomials.VanderDerivatives(reference_element_data.Monomials,
                                                                           localSpace.VanderBoundary,
                                                                           polygonDiameter);

        // Compute mass matrix of monomials.
        localSpace.Hmatrix = localSpace.VanderInternal.transpose() *
                             internalQuadratureWeights.asDiagonal() *
                             localSpace.VanderInternal;
        // Compute LLT factorization of order-1 monomials.
        localSpace.H_km1_LLT = localSpace.Hmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).llt();

        localSpace.Qmatrix = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions,
                                                localSpace.NumProjectorBasisFunctions);
        localSpace.QmatrixInv = MatrixXd::Identity(localSpace.NumProjectorBasisFunctions,
                                                   localSpace.NumProjectorBasisFunctions);
      }
      //****************************************************************************
      Eigen::MatrixXd VEM_PCC_2D_LocalSpace::ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                            const double& polygonMeasure,
                                                            const double& polygonDiameter,
                                                            const Eigen::VectorXd& internalQuadratureWeights,
                                                            const Eigen::VectorXd& boundaryQuadratureWeights,
                                                            const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                            VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        // G_{ij} = \int_E \nabla m_i \nabla m_j
        MatrixXd Gmatrix = localSpace.VanderInternalDerivatives[0].transpose() *
                           internalQuadratureWeights.asDiagonal() *
                           localSpace.VanderInternalDerivatives[0] +
                           localSpace.VanderInternalDerivatives[1].transpose() *
                           internalQuadratureWeights.asDiagonal() *
                           localSpace.VanderInternalDerivatives[1];
        // B_{ij} = \int_E \nabla m_i \nabla \phi_j
        MatrixXd Bmatrix;
        Bmatrix.setZero(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);
        // First block of B: \int_{\partial E}\frac{\partial m_i}{\partial n} \phi_j
        Bmatrix.leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundaryDerivatives[0].transpose() *
            boundaryQuadratureWeightsTimesNormal[0].asDiagonal() +
            localSpace.VanderBoundaryDerivatives[1].transpose()*
            boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

        if (localSpace.Order == 1)
        {
          // B_{0j} = \int_{\partial E} \phi_j
          Bmatrix.row(0) = boundaryQuadratureWeights;
          // G_{0j} = \int_{\partial E} m_j
          Gmatrix.row(0) = localSpace.VanderBoundary.transpose() *
                           boundaryQuadratureWeights;
        }
        else
        {
          // G_{0j} = \int_{E} m_j
          Gmatrix.row(0) = localSpace.VanderInternal.transpose() *
                           internalQuadratureWeights;
          // Second block of B: - \int_E \Delta m_i \phi_j
          Bmatrix.rightCols(localSpace.NumInternalBasisFunctions) =
              (- polygonMeasure / (polygonDiameter * polygonDiameter)) *
              (reference_element_data.Monomials.Laplacian.leftCols(localSpace.NumInternalBasisFunctions));
          // B_{0j} = \int_{E} \phi_j (only the first internal basis
          // function has a non-zero integral)
          Bmatrix(0, localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) = polygonMeasure;
        }

        localSpace.Bmatrix = Bmatrix;
        localSpace.Gmatrix = Gmatrix;

        return Gmatrix.partialPivLu().solve(Bmatrix);
      }
      //****************************************************************************
      MatrixXd VEM_PCC_2D_LocalSpace::ComputeBasisPolynomialsDofs(const double& polygonMeasure,
                                                                  const VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        MatrixXd basisPolynomialsDofs;
        basisPolynomialsDofs.resize(localSpace.NumBasisFunctions, localSpace.NumProjectorBasisFunctions);
        // boundary degrees of freedom of monomials (values at points on
        // the boundary)
        basisPolynomialsDofs.topRows(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary;

        if (localSpace.Order > 1)
        {
          // internal degrees of freedom of monomials (scaled moments)
          basisPolynomialsDofs.bottomRows(localSpace.NumInternalBasisFunctions) =
              localSpace.Hmatrix.topRows(localSpace.NumInternalBasisFunctions) / polygonMeasure;
        }

        return basisPolynomialsDofs;
      }
      //****************************************************************************
      void VEM_PCC_2D_LocalSpace::ComputeL2Projectors(const double& polygonMeasure,
                                                      VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        Eigen::MatrixXd Cmatrix(localSpace.NumProjectorBasisFunctions, localSpace.NumBasisFunctions);
        // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
        Cmatrix.bottomRows(localSpace.NumProjectorBasisFunctions - localSpace.NumInternalBasisFunctions) =
            localSpace.Hmatrix.bottomRows(localSpace.NumProjectorBasisFunctions - localSpace.NumInternalBasisFunctions) *
            localSpace.PiNabla;

        if (localSpace.Order > 1)
        {
          Cmatrix.topLeftCorner(localSpace.NumInternalBasisFunctions,
                                localSpace.NumBasisFunctions - localSpace.NumInternalBasisFunctions).setZero();
          //\int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
          Cmatrix.topRightCorner(localSpace.NumInternalBasisFunctions, localSpace.NumInternalBasisFunctions) =
              polygonMeasure * Eigen::MatrixXd::Identity(localSpace.NumInternalBasisFunctions,
                                                         localSpace.NumInternalBasisFunctions);
        }

        localSpace.Pi0km1 = localSpace.H_km1_LLT.solve(Cmatrix.topRows(localSpace.Nkm1));
        localSpace.Pi0k = localSpace.Hmatrix.llt().solve(Cmatrix);
        localSpace.Cmatrix = Cmatrix;
      }
      //****************************************************************************
      void VEM_PCC_2D_LocalSpace::ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                   const double& polygonMeasure,
                                                                   const double& polygonDiameter,
                                                                   const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                                   VEM_PCC_2D_LocalSpace_Data& localSpace) const
      {
        localSpace.Ematrix.resize(2, MatrixXd::Zero(localSpace.Nkm1, localSpace.NumBasisFunctions));
        localSpace.Ematrix[0].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() * boundaryQuadratureWeightsTimesNormal[0].asDiagonal();

        localSpace.Ematrix[1].leftCols(localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions) =
            localSpace.VanderBoundary.leftCols(localSpace.Nkm1).transpose() *
            boundaryQuadratureWeightsTimesNormal[1].asDiagonal();

        if (localSpace.Order > 1)
        {
          localSpace.Ematrix[0].rightCols(localSpace.NumInternalBasisFunctions) =
              -(polygonMeasure / polygonDiameter) *
              monomials.D_x(reference_element_data.Monomials).topLeftCorner(localSpace.Nkm1,
                                                                            localSpace.NumInternalBasisFunctions);
          localSpace.Ematrix[1].rightCols(localSpace.NumInternalBasisFunctions) =
              -(polygonMeasure / polygonDiameter) *
              monomials.D_y(reference_element_data.Monomials).topLeftCorner(localSpace.Nkm1,
                                                                            localSpace.NumInternalBasisFunctions);
        }

        localSpace.Pi0km1Der.resize(2);
        localSpace.Pi0km1Der[0] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[0]);
        localSpace.Pi0km1Der[1] = localSpace.H_km1_LLT.solve(localSpace.Ematrix[1]);
      }
    }
  }
}
