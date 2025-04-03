#ifndef __I_VEM_MCC_2D_Pressure_LocalSpace_HPP
#define __I_VEM_MCC_2D_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
/// \brief Interface class for the pressure space of Mixed Conforming Constant degree 2D \cite secondMixed \cite
/// DaVeiga2016 \cite Teora2023 \cite Teora2024_mixed.
class I_VEM_MCC_2D_Pressure_LocalSpace
{
  public:
    /// \brief Class destructor.
    virtual ~I_VEM_MCC_2D_Pressure_LocalSpace()
    {
    }

    /// \brief Compute data of VEM space on a polygon.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \return A data structure containing the information about the local space.
    virtual VEM_MCC_2D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_MCC_2D_Polygon_Geometry &polygon) const = 0;

    /// \brief Compute values of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of basis functions at given points.
    /// \param reference_element_data Data of the reference element.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                        const VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                        const Eigen::MatrixXd &points) const = 0;
    /// \brief Compute values of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                     const VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
