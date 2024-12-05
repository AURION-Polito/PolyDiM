#ifndef __VEM_MCC_2D_VelocityLocalSpace_HPP
#define __VEM_MCC_2D_VelocityLocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_MCC_VelocityLocalSpace_Data.hpp"
#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
/// \brief Class used for computing values of basis functions of 2D
/// Mixed Conforming Constant degree Virtual Element Methods.
class VEM_MCC_2D_VelocityLocalSpace final
{
private:
    MCC::VEM_MCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                         const Eigen::MatrixXd& polygonVertices,
                                         const Eigen::Vector3d& polygonCentroid,
                                         const double& polygonDiameter,
                                         const Eigen::MatrixXd& internalQuadraturePoints,
                                         const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::MatrixXd& boundaryQuadraturePoints,
                                         VEM_MCC_VelocityLocalSpace_Data& localSpace) const;

    void ComputeL2Projectors(const double& polygonMeasure,
                             VEM_MCC_VelocityLocalSpace_Data& localSpace) const;


    inline Eigen::MatrixXd ComputeStabilizationMatrixPi0k(const double& polygonMeasure,
                                                          const VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        return utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                        polygonMeasure,
                                                        localSpace.Dmatrix);
    }

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       VEM_MCC_VelocityLocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const Eigen::MatrixXd &polytopeVertices,
                                 const Eigen::MatrixXd &edgeNormals,
                                 const std::vector<bool> &edgeDirections,
                                 const Eigen::VectorXd &boundaryQuadratureWeights,
                                 Eigen::VectorXd &edgeDirectionsVector,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 VEM_MCC_VelocityLocalSpace_Data &localSpace) const;
public:
    VEM_MCC_VelocityLocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                        const VEM_MCC_2D_Polygon_Geometry& polygon) const;

    /// \brief Compute matrix D: D_{ij} = dof_i(m_j).
    Eigen::MatrixXd ComputePolynomialBasisDofs(const double& polytopeMeasure,
                                               const VEM_MCC_VelocityLocalSpace_Data& localSpace) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_VelocityLocalSpace_Data& localSpace,
                                                       const ProjectionTypes& projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_MCC_VelocityLocalSpace_Data& localSpace,
                                                                              const ProjectionTypes& projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.VanderInternalDerivatives,
                                                              localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                       const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                       const VEM_MCC_VelocityLocalSpace_Data& localSpace,
                                                       const ProjectionTypes& projectionType,
                                                       const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data,
                                                                              polygon,
                                                                              points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                              const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                                              const VEM_MCC_VelocityLocalSpace_Data& localSpace,
                                                                              const ProjectionTypes& projectionType,
                                                                              const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               ComputePolynomialsValues(reference_element_data,
                                                                                        polygon,
                                                                                        points),
                                                               ComputePolynomialsDerivativeValues(reference_element_data,
                                                                                                  polygon,
                                                                                                  points),
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                                const VEM_MCC_VelocityLocalSpace_Data& localSpace,
                                                                const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.Pi0km1Der,
                                                              ComputePolynomialsDerivativeValues(reference_element_data,
                                                                                                 polygon,
                                                                                                 points));
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                    const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                    const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  polygon.Centroid,
                                                  polygon.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                           const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                                           const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polygon.Diameter,
                                                            ComputePolynomialsValues(reference_element_data,
                                                                                     polygon,
                                                                                     points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                             const VEM_MCC_2D_Polygon_Geometry& polygon,
                                                             const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Diameter,
                                                           ComputePolynomialsValues(reference_element_data,
                                                                                    polygon,
                                                                                    points));
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_MCC_VelocityLocalSpace_Data& localSpace,
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
