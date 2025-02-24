#ifndef __I_VEM_DF_PCC_2D_Velocity_LocalSpace_HPP
#define __I_VEM_DF_PCC_2D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Class used for computing values of basis functions of 2D
/// Divergence Free Primal Conforming Constant degree Virtual Element Methods.
class I_VEM_DF_PCC_2D_Velocity_LocalSpace
{
  public:
    virtual ~I_VEM_DF_PCC_2D_Velocity_LocalSpace()
    {
    }

    virtual VEM_DF_PCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const = 0;

    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeValuesOnEdge(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                     const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                     const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType,
                                                                     const Eigen::MatrixXd &points) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                               const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                               const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType,
                                                                               const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                     const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                     const Eigen::MatrixXd &points) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                            const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                              const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                              const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
