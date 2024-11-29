#ifndef __VEM_PCC_2D_LocalSpace_HPP
#define __VEM_PCC_2D_LocalSpace_HPP

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
      /// \brief Class used for computing values of basis functions of 2D
      /// Primal Conforming Constant degree Virtual Element Methods.
      class VEM_PCC_2D_LocalSpace final
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
          /// \note This requires \ref InitializeProjectorsComputation() to be called previously.
          Eigen::MatrixXd ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                         const double& polygonMeasure,
                                         const double& polygonDiameter,const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::VectorXd& boundaryQuadratureWeights,
                                         const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                         VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrices \ref pi0km1 and \ref pi0k.
          /// \note This requires \ref ComputePiNabla() to be called previously.
          void ComputeL2Projectors(const double& polygonMeasure,
                                   VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          /// \brief Compute matrices \ref pi0km1Der.
          void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                const double& polygonMeasure,
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
                                                        ComputeBasisPolynomialsDofs(polygonMeasure,
                                                                                    localSpace));
          }
          /// \brief Compute matrix \ref stabMatrix with Pi0k projector.
          /// \note This requires \ref ComputeL2Projectors() to be called previously.
          inline Eigen::MatrixXd ComputeStabilizationMatrixPi0k(const double& polygonMeasure,
                                                                const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                            polygonMeasure,
                                                            ComputeBasisPolynomialsDofs(polygonMeasure,
                                                                                        localSpace));
          }

        public:
          VEM_PCC_2D_LocalSpace();
          virtual ~VEM_PCC_2D_LocalSpace() {}

          VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                      const VEM_PCC_2D_Polygon_Geometry& polygon) const;

          VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                        const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                        const Eigen::MatrixXd& internalQuadraturePoints,
                                                        const Eigen::VectorXd& internalQuadratureWeights,
                                                        const Eigen::MatrixXd& boundaryQuadraturePoints,
                                                        const Eigen::VectorXd& boundaryQuadratureWeights,
                                                        const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal) const;

          /// \brief Compute matrix D: D_{ij} = dof_i(m_j).
          Eigen::MatrixXd ComputeBasisPolynomialsDofs(const double& polytopeMeasure,
                                                      const VEM_PCC_2D_LocalSpace_Data& localSpace) const;

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                             const ProjectionTypes& projectionType) const
          {
            return utilities.ComputeBasisFunctionsValues(projectionType,
                                                         localSpace.Nkm1,
                                                         localSpace.Pi0km1,
                                                         localSpace.Pi0k,
                                                         localSpace.VanderInternal);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                                    const ProjectionTypes& projectionType) const
          {
            return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                                   localSpace.Nkm1,
                                                                   localSpace.VanderInternal,
                                                                   localSpace.VanderInternalDerivatives,
                                                                   localSpace.PiNabla,
                                                                   localSpace.Pi0km1Der);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                                  localSpace.VanderInternalDerivatives,
                                                                  localSpace.Pi0km1Der);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                             const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                             const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                             const ProjectionTypes& projectionType,
                                                             const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisFunctionsValues(projectionType,
                                                         localSpace.Nkm1,
                                                         localSpace.Pi0km1,
                                                         localSpace.Pi0k,
                                                         ComputeBasisPolynomialsValues(reference_element_data,
                                                                                       polygon,
                                                                                       points));
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                                    const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                                                    const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                                    const ProjectionTypes& projectionType,
                                                                                    const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                                   localSpace.Nkm1,
                                                                   ComputeBasisPolynomialsValues(reference_element_data,
                                                                                                 polygon,
                                                                                                 points),
                                                                   ComputeBasisPolynomialsDerivativeValues(reference_element_data,
                                                                                                           polygon,
                                                                                                           points),
                                                                   localSpace.PiNabla,
                                                                   localSpace.Pi0km1Der);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                      const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                                      const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                                      const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                                  localSpace.Pi0km1Der,
                                                                  ComputeBasisPolynomialsDerivativeValues(reference_element_data,
                                                                                                          polygon,
                                                                                                          points));
          }

          inline Eigen::MatrixXd ComputeBasisPolynomialsValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisPolynomialsValues(localSpace.VanderInternal);
          }

          inline Eigen::MatrixXd ComputeBasisPolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                               const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                               const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisPolynomialsValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Centroid,
                                                           polygon.Diameter,
                                                           points);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisPolynomialsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data& localSpace) const
          {
            return utilities.ComputeBasisPolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisPolynomialsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                                      const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                                                      const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisPolynomialsDerivativeValues(reference_element_data.Monomials,
                                                                     monomials,
                                                                     polygon.Diameter,
                                                                     ComputeBasisPolynomialsValues(reference_element_data,
                                                                                                   polygon,
                                                                                                   points));
          }

          inline Eigen::MatrixXd ComputeBasisPolynomialsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                                        const VEM_PCC_2D_Polygon_Geometry& polygon,
                                                                        const Eigen::MatrixXd& points) const
          {
            return utilities.ComputeBasisPolynomialsLaplacianValues(reference_element_data.Monomials,
                                                                    monomials,
                                                                    polygon.Diameter,
                                                                    ComputeBasisPolynomialsValues(reference_element_data,
                                                                                                  polygon,
                                                                                                  points));
          }

          inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_LocalSpace_Data& localSpace,
                                                     const Eigen::VectorXd& edgeInternalPoints,
                                                     const Eigen::VectorXd& pointsCurvilinearCoordinates) const
          {
            const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(localSpace.Order,
                                                                                                 edgeInternalPoints);
            return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                                 localSpace.Order,
                                                 edgeBasisCoefficients,
                                                 pointsCurvilinearCoordinates);
          }
      };
    }
  }
}

#endif
