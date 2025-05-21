// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

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
/// \brief Interface class for the pressure local space for 3D Mixed Conforming Constant degree Virtual Element Methods.
class I_VEM_MCC_3D_Pressure_LocalSpace
{
  public:
    /// \brief Class destructor.
    virtual ~I_VEM_MCC_3D_Pressure_LocalSpace()
    {
    }

    /// \brief Compute data of pressure space on a polygon.
    /// \param reference_element_data Data of the reference element of the pressure space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \return A data structure containing the information about the local pressure space.
    virtual VEM_MCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_MCC_3D_Polyhedron_Geometry &polyhedron) const = 0;

    /// \brief Compute values of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix with values of a basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of basis functions at given points.
    /// \param reference_element_data Data of the reference element of the pressure space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix with values of a basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                        const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                     const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
