#ifndef __VEM_PCC_2D_LocalSpace_HPP
#define __VEM_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim {
namespace VEM {
namespace PCC {
/// \brief Class used for computing values of basis functions of 2D
/// Primal Conforming Constant degree Virtual Element Methods.
class VEM_PCC_2D_LocalSpace final
{
private:
    VEM_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const Eigen::VectorXd &boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                        VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                          const double &polygonMeasure,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polygonMeasure,
                             VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        utilities.ComputeL2Projectors(polygonMeasure,
                                      localSpace.Order,
                                      localSpace.Nkm1,
                                      localSpace.NumProjectorBasisFunctions,
                                      localSpace.NumInternalBasisFunctions,
                                      localSpace.NumBasisFunctions,
                                      localSpace.Hmatrix,
                                      localSpace.PiNabla,
                                      localSpace.Cmatrix,
                                      localSpace.Pi0km1,
                                      localSpace.Pi0k);
    };

    void ComputePolynomialsDofs(const double &polytopeMeasure,
                                VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    inline void ComputeStabilizationMatrix(const double &polygonDiameter,
                                           VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        localSpace.StabMatrix = utilities.ComputeStabilizationMatrix(localSpace.PiNabla,
                                                                     polygonDiameter,
                                                                     localSpace.Dmatrix);
    }

    inline void ComputeStabilizationMatrixPi0k(const double &polygonMeasure,
                                               VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        localSpace.StabMatrixPi0k = utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                                             polygonMeasure,
                                                                             localSpace.Dmatrix);
    }

public:

    /// \brief Create and Initialize all the variables contained in \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data
    /// \input The functions requires
    ///         - an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains monomials, quadrature and the number of degrees of freedom,
    /// counting in order DOFS associated with vertices, edges and internal values.
    ///         - an object of type \ref VEM::PCC::VEM_PCC_2D_Polygon_Geometry which contains the geoemtric properties of the elements.
    /// \return An object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data.
    VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                  const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                  const Eigen::MatrixXd &internalQuadraturePoints,
                                                  const Eigen::VectorXd &internalQuadratureWeights,
                                                  const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                  const Eigen::VectorXd &boundaryQuadratureWeights,
                                                  const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.VanderInternalDerivatives,
                                                              localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                       const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data,
                                                                              polygon,
                                                                              points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                                              const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, polygon, points),
            ComputePolynomialsDerivativeValues(reference_element_data, polygon, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                                const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.Pi0km1Der,
                                                              ComputePolynomialsDerivativeValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  polygon.Centroid,
                                                  polygon.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polygon.Diameter,
                                                            ComputePolynomialsValues(reference_element_data,
                                                                                     polygon,
                                                                                     points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const VEM_PCC_2D_Polygon_Geometry &polygon,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Diameter,
                                                           ComputePolynomialsValues(reference_element_data,
                                                                                    polygon,
                                                                                    points));
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge( const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                               const Eigen::VectorXd &edgeInternalPoints,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(localSpace.Order, edgeInternalPoints);
        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                             localSpace.Order,
                                             edgeBasisCoefficients,
                                             pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
