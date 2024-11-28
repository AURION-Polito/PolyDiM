#ifndef __VEM_PCC_2D_LocalSpace_HPP
#define __VEM_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
  namespace VEM
  {
    namespace PCC
    {
      struct VEM_PCC_2D_Polygon_Geometry final
      {
          Eigen::MatrixXd& Vertices;
          Eigen::Vector3d& Centroid;
          double& Measure;
          double& Diameter;
      };

      /// \brief Class used for computing values of basis functions of 2D
      /// Primal Conforming Constant degree Virtual Element Methods.
      class VEM_PCC_2D_LocalSpace final
      {
        private:
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
          void InitializeProjectorsComputation(const Eigen::MatrixXd& polygonVertices,
                                               const Eigen::Vector3d& polygonCentroid,
                                               const double& polygonDiameter,
                                               const Eigen::MatrixXd& internalQuadraturePoints,
                                               const Eigen::VectorXd& internalQuadratureWeights,
                                               const Eigen::MatrixXd& boundaryQuadraturePoints,
                                               VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrix \ref piNabla.
          /// \note This requires \ref InitializeProjectorsComputation() to be called previously.
          Eigen::MatrixXd ComputePiNabla(const double& polygonMeasure,
                                         const double& polygonDiameter,const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::VectorXd& boundaryQuadratureWeights,
                                         const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                         VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrices \ref pi0km1 and \ref pi0k.
          /// \note This requires \ref ComputePiNabla() to be called previously.
          void ComputeL2Projectors(const double& polygonMeasure,
                                   VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrices \ref pi0km1Der.
          void ComputeL2ProjectorsOfDerivatives(const double& polygonMeasure,
                                                const double& polygonDiameter,
                                                const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                                VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute the stabilization matrix with PiNabla projector.
          /// \note used for method with stabilization
          /// \return stabilization matrix, size numQuadraturePoints x NumberBasisFunctions()
          inline Eigen::MatrixXd ComputeStabilizationMatrix(const double& polygonMeasure,
                                                            const double& polygonDiameter,
                                                            VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeStabilizationMatrix(localSpace.PiNabla,
                                                        polygonDiameter,
                                                        ComputePolynomialBasisDofs(polygonMeasure,
                                                                                   localSpace));
          }
          /// \brief Compute matrix \ref stabMatrix with Pi0k projector.
          /// \note This requires \ref ComputeL2Projectors() to be called previously.
          inline Eigen::MatrixXd ComputeStabilizationMatrixPi0k(const double& polygonMeasure,
                                                                const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                            polygonMeasure,
                                                            ComputePolynomialBasisDofs(polygonMeasure,
                                                                                       localSpace));
          }

        public:



          VEM_PCC_2D_LocalSpace();
          virtual ~VEM_PCC_2D_LocalSpace() {}

          VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const unsigned int order,
                                                      const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                      const Eigen::MatrixXd& internalQuadraturePoints,
                                                      const Eigen::VectorXd& internalQuadratureWeights,
                                                      const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                      const Eigen::VectorXd& boundaryQuadratureWeights,
                                                      const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal) const;

          VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                        const Eigen::MatrixXd& internalQuadraturePoints,
                                                        const Eigen::VectorXd& internalQuadratureWeights,
                                                        const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                        const Eigen::VectorXd& boundaryQuadratureWeights,
                                                        const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal) const;

          /// \brief Compute matrix D: D_{ij} = dof_i(m_j).
          Eigen::MatrixXd ComputePolynomialBasisDofs(const double& polytopeMeasure,
                                                     const VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          inline Eigen::MatrixXd ComputeBasisFunctionValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisFunctionValues(localSpace.ProjectionType == VEM_PCC_2D_LocalSpace_Data::ProjectionTypes::Pi0km1,
                                                        localSpace.Nkm1,
                                                        localSpace.Pi0km1,
                                                        localSpace.Pi0k,
                                                        localSpace.VanderInternal);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisFunctionDerivativeValues(localSpace.Nkm1,
                                                                  localSpace.VanderInternal,
                                                                  localSpace.Pi0km1Der);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValuesPiNabla(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues;

            basisFunctionDerivativeValues.resize(referenceElement.Dimension());
            for(unsigned short i = 0; i < referenceElement.Dimension(); ++i)
              basisFunctionDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
          }

          inline Eigen::MatrixXd ComputeBasisFunctionLaplacianValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisFunctionLaplacianValues(localSpace.Nkm1,
                                                                 localSpace.VanderInternalDerivatives,
                                                                 localSpace.Pi0km1Der);
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
            return utilities.ComputeBasisFunctionLaplacianValuesOnPoints(localSpace.Nkm1,
                                                                         localSpace.Pi0km1Der,
                                                                         ComputePolynomialBasisDerivativeValues(polytopeCentroid,
                                                                                                                polytopeDiameter,
                                                                                                                localSpace,
                                                                                                                points));
          }

          inline Eigen::MatrixXd ComputePolynomialBasisValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputePolynomialBasisValues(localSpace.VanderInternal);
          }

          inline Eigen::MatrixXd ComputePolynomialBasisValues(const Eigen::Vector3d& polytopeCentroid,
                                                              const double& polytopeDiameter,
                                                              const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                              const Eigen::MatrixXd& points) const
          {
            return utilities.ComputePolynomialBasisValues(monomials,
                                                          polytopeCentroid,
                                                          polytopeDiameter,
                                                          points);
          }

          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputePolynomialBasisDerivativeValues(localSpace.VanderInternalDerivatives);
          }

          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const Eigen::Vector3d& polytopeCentroid,
                                                                                     const double& polytopeDiameter,
                                                                                     const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                                     const Eigen::MatrixXd& points) const
          {
            return utilities.ComputePolynomialBasisDerivativeValues(monomials,
                                                                    polytopeDiameter,
                                                                    ComputePolynomialBasisValues(polytopeCentroid,
                                                                                                 polytopeDiameter,
                                                                                                 localSpace,
                                                                                                 points));
          }

          inline Eigen::MatrixXd ComputePolynomialBasisLaplacianValues(const Eigen::Vector3d& polytopeCentroid,
                                                                       const double& polytopeDiameter,
                                                                       const Eigen::MatrixXd& internalQuadraturePoints,
                                                                       const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return ComputePolynomialBasisLaplacianValues(polytopeCentroid,
                                                         polytopeDiameter,
                                                         localSpace,
                                                         internalQuadraturePoints);
          }

          inline Eigen::MatrixXd ComputePolynomialBasisLaplacianValues(const Eigen::Vector3d& polytopeCentroid,
                                                                       const double& polytopeDiameter,
                                                                       const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                       const Eigen::MatrixXd& points) const
          {
            return utilities.ComputePolynomialBasisLaplacianValues(monomials,
                                                                   polytopeDiameter,
                                                                   ComputePolynomialBasisValues(polytopeCentroid,
                                                                                                polytopeDiameter,
                                                                                                localSpace,
                                                                                                points));
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
