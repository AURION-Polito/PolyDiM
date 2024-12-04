#ifndef __VEM_PCC_2D_Ortho_LocalSpace_HPP
#define __VEM_PCC_2D_Ortho_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
  namespace VEM
  {
    namespace PCC
    {
      class VEM_PCC_2D_Ortho_LocalSpace final
      {
        private:
          VEM_PCC_Utilities<2> utilities;
          Monomials::VEM_Monomials_2D monomials;

          /// \brief Initialize quantities required for computing projectors.
          /// \details This method computes \ref measure, \ref diameter, \ref
          /// vanderInternal, \ref vanderInternalDerivatives, \ref
          /// internalWeights, \ref vanderBoundary, \ref
          /// vanderBoundaryDerivatives, \ref boundaryWeights, \ref
          /// boundaryWeightsTimesNormal, \ref Hmatrix.
          /// \param geometry The geometry used as domain for the computation of projectors. It has to be
          /// an object of class \ref Polygon with dimension 2.
          /// \note The following methods have to be called on the geometry before calling this method:
          ///  - \ref Polygon::Compute2DPolygonProperties()
          ///  - \ref Polygon::ComputePositionPoint()
          ///  - \ref Polygon::ComputeNormalSign()
          ///  .
          /// Moreover, the method \ref Segment::ComputeNormal() has to be called on each edge of the
          /// geometry.
          void InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                               const Eigen::MatrixXd& polygonVertices,
                                               const Eigen::Vector3d& polygonCentroid,
                                               const double& polygonDiameter,
                                               const Eigen::MatrixXd& internalQuadraturePoints,
                                               const Eigen::VectorXd& internalQuadratureWeights,
                                               const Eigen::MatrixXd& boundaryQuadraturePoints,
                                               VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrix \ref piNabla.
          /// \note This requires \ref InitializeProjectorsComputation() to be
          /// called previously.
          Eigen::MatrixXd ComputePiNabla(const double& polygonMeasure,
                                         const double& polygonDiameter,
                                         const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::VectorXd& boundaryQuadratureWeights,
                                         const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                         VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          void ComputeL2Projectors(const double& polygonMeasure,
                                   VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          inline Eigen::MatrixXd ComputeStabilizationMatrix(const double& polygonMeasure,
                                                            const double& polygonDiameter,
                                                            VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeStabilizationMatrix(localSpace.PiNabla,
                                                        polygonDiameter,
                                                        ComputePolynomialsDofs(polygonMeasure,
                                                                               localSpace));
          }

          /// \brief Compute matrix \ref stabMatrix with Pi0k projector.
          /// \note This requires \ref ComputeL2Projectors() to be called previously.
          inline Eigen::MatrixXd ComputeStabilizationMatrixPi0k(const double& polygonMeasure,
                                                                const double& polygonDiameter,
                                                                const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                            polygonDiameter,
                                                            ComputePolynomialsDofs(polygonMeasure,
                                                                                   localSpace));
          }

          /// Compute the change of basis matrix and the mass matrix of orthogonal polynomial basis
          void ChangeOfBasis(const Eigen::VectorXd& internalQuadratureWeights,
                             VEM_PCC_2D_LocalSpace_Data& localSpace) const;

        public:
          VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                      const VEM_PCC_2D_Polygon_Geometry& polygon) const;

          VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const Eigen::MatrixXd& polygonVertices,
                                                        const Eigen::Vector3d& polygonCentroid,
                                                        const double& polygonMeasure,
                                                        const double& polygonDiameter,
                                                        const Eigen::MatrixXd& internalQuadraturePoints,
                                                        const Eigen::VectorXd& internalQuadratureWeights,
                                                        const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                        const Eigen::VectorXd& boundaryQuadratureWeights,
                                                        const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                        const VEM_PCC_2D_LocalSpace_Data::ProjectionTypes& projectionType) const;

          /// \brief Compute matrices \ref pi0km1Der.
          /// \note This requires \ref InitializeProjectorsComputation() to
          /// be called previously..
          void ComputeL2ProjectorsOfDerivatives(const double& polygonMeasure,
                                                const double& polygonDiameter,
                                                const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          Eigen::MatrixXd ComputePolynomialsDofs(const double& polytopeMeasure,
                                                 const VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          inline Eigen::MatrixXd ComputeBasisFunctionValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            switch (localSpace.ProjectionType)
            {
              case VEM_PCC_2D_LocalSpace_Data::ProjectionTypes::Pi0km1:
                return localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                    localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                     localSpace.Nkm1).transpose() *
                    localSpace.Pi0km1;
              case VEM_PCC_2D_LocalSpace_Data::ProjectionTypes::Pi0k:
                return localSpace.VanderInternal *
                    localSpace.Qmatrix.transpose() *
                    localSpace.Pi0k;
              default:
                throw std::runtime_error("Unsupported ProjectionTypes");
            }
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues;

            basisFunctionsDerivativeValues.resize(referenceElement.Dimension());
            for(unsigned short i = 0; i < referenceElement.Dimension(); ++i)
              basisFunctionsDerivativeValues[i] = localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                                                  localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                   localSpace.Nkm1).transpose() *
                                                  localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValuesPiNabla(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues;

            basisFunctionDerivativeValues.resize(referenceElement.Dimension());
            for(unsigned short i = 0; i < referenceElement.Dimension(); ++i)
              basisFunctionDerivativeValues[i] = localSpace.VanderInternalDerivatives[i]
                                                 * localSpace.Qmatrix.transpose()
                                                 * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
          }

          inline Eigen::MatrixXd ComputeBasisFunctionLaplacianValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            throw std::runtime_error("Unimplemented method");
          }

          inline Eigen::MatrixXd ComputeBasisFunctionValues(const Eigen::Vector3d& polytopeCentroid,
                                                            const double& polytopeDiameter,
                                                            const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                            const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisFunctionValues(localSpace.ProjectionType == VEM_PCC_2D_LocalSpace_Data::ProjectionTypes::Pi0km1,
                                                        localSpace.Nkm1,
                                                        localSpace.Pi0km1,
                                                        localSpace.Pi0k,
                                                        localSpace.VanderInternal,
                                                        ComputePolynomialBasisValues(polytopeCentroid,
                                                                                     polytopeDiameter,
                                                                                     localSpace,
                                                                                     points));
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const Eigen::Vector3d& polytopeCentroid,
                                                                                   const double& polytopeDiameter,
                                                                                   const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                                   const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisFunctionDerivativeValues(localSpace.Nkm1,
                                                                  localSpace.Pi0km1Der,
                                                                  ComputePolynomialBasisValues(polytopeCentroid,
                                                                                               polytopeDiameter,
                                                                                               localSpace,
                                                                                               points));
          }

          inline Eigen::MatrixXd ComputeBasisFunctionLaplacianValues(const Eigen::Vector3d& polytopeCentroid,
                                                                     const double& polytopeDiameter,
                                                                     const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                     const Eigen::MatrixXd& points) const
          {
            throw std::runtime_error("Unimplemented method");
          }

          inline Eigen::MatrixXd ComputePolynomialBasisValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
          }

          inline Eigen::MatrixXd ComputePolynomialBasisValues(const Eigen::Vector3d& polytopeCentroid,
                                                              const double& polytopeDiameter,
                                                              const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                              const Eigen::MatrixXd& points) const
          {
            return monomials.Vander(points,
                                    polytopeCentroid,
                                    polytopeDiameter) * localSpace.Qmatrix.transpose();
          }

          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            vector<Eigen::MatrixXd> polynomialBasisDerivativeValues;
            polynomialBasisDerivativeValues.resize(referenceElement.Dimension());
            for(unsigned short i = 0; i < referenceElement.Dimension(); ++i)
              polynomialBasisDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] *
                                                   localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                    localSpace.Nkm1).transpose();
            return polynomialBasisDerivativeValues;
          }

          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const Eigen::Vector3d& polytopeCentroid,
                                                                                     const double& polytopeDiameter,
                                                                                     const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                                     const Eigen::MatrixXd& points) const
          {
            Eigen::MatrixXd monomialBasisValues;
            vector<Eigen::MatrixXd> monomialBasisDerivativeValues;
            monomialBasisValues = monomials.Vander(points,
                                                   polytopeCentroid,
                                                   polytopeDiameter);
            monomialBasisDerivativeValues = monomials.VanderDerivatives(monomialBasisValues,
                                                                        polytopeDiameter);

            vector<Eigen::MatrixXd> polynomialBasisDerivativeValues;
            polynomialBasisDerivativeValues.resize(referenceElement.Dimension());
            for(unsigned short i = 0; i < referenceElement.Dimension(); ++i)
              polynomialBasisDerivativeValues[i] = monomialBasisDerivativeValues[i] *
                                                   localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                    localSpace.Nkm1).transpose();
            return polynomialBasisDerivativeValues;
          }

          inline Eigen::MatrixXd ComputePolynomialBasisLaplacianValues(const Eigen::Vector3d& polytopeCentroid,
                                                                       const double& polytopeDiameter,
                                                                       const Eigen::MatrixXd& internalQuadraturePoints,
                                                                       const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            throw std::runtime_error("Unimplemented method");
          }

          inline Eigen::MatrixXd ComputePolynomialBasisLaplacianValues(const Eigen::Vector3d& polytopeCentroid,
                                                                       const double& polytopeDiameter,
                                                                       const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                       const Eigen::MatrixXd& points) const
          {
            throw std::runtime_error("Unimplemented method");
          }

          inline Eigen::MatrixXd ComputeValuesOnEdge(const Eigen::VectorXd& edgeInternalPoints,
                                                     const Eigen::VectorXd& pointsCurvilinearCoordinates) const
          {
            const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(referenceElement.Order(),
                                                                                                 edgeInternalPoints);
            return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                                 referenceElement.Order(),
                                                 edgeBasisCoefficients,
                                                 pointsCurvilinearCoordinates);
          }
      };
    }
  }
}

#endif
