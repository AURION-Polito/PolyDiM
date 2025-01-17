#ifndef __I_VEM_MCC_2D_Velocity_LocalSpace_HPP
#define __I_VEM_MCC_2D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
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
class I_VEM_MCC_2D_Velocity_LocalSpace
{
  public:
    virtual ~I_VEM_MCC_2D_Velocity_LocalSpace()
    {
    }

    virtual VEM_MCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_MCC_2D_Polygon_Geometry &polygon) const = 0;

    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                     const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
