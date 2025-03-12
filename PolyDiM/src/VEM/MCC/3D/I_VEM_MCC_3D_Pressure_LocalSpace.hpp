#ifndef __I_VEM_MCC_3D_Pressure_LocalSpace_HPP
#define __I_VEM_MCC_3D_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_MCC_3D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{

class I_VEM_MCC_3D_Pressure_LocalSpace
{
  public:
  /// \brief Class destructor.
  virtual ~I_VEM_MCC_3D_Pressure_LocalSpace()
    {
    }

    virtual VEM_MCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_MCC_3D_Polyhedron_Geometry &polyhedron) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                        const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                        const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                     const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
