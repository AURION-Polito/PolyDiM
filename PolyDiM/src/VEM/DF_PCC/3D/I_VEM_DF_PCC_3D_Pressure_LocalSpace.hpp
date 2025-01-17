#ifndef __I_VEM_DF_PCC_3D_Pressure_LocalSpace_HPP
#define __I_VEM_DF_PCC_3D_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class I_VEM_DF_PCC_3D_Pressure_LocalSpace
{
  public:
    virtual ~I_VEM_DF_PCC_3D_Pressure_LocalSpace()
    {
    }

    virtual VEM_DF_PCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                        const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                        const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                     const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
