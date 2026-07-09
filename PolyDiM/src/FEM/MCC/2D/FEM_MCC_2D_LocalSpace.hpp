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

#ifndef __FEM_MCC_2D_LocalSpace_HPP
#define __FEM_MCC_2D_LocalSpace_HPP

#include "FEM_MCC_2D_LocalSpace_Data.hpp"
#include "FEM_MCC_2D_ReferenceElement.hpp"
#include "FEM_Triangle_RT_MCC_2D_LocalSpace.hpp"

namespace Polydim
{
namespace FEM
{
namespace MCC
{

/// @brief Local space for a 2D mixed (velocity/pressure) MCC finite element.
///
/// Provides a uniform, type-agnostic interface to build the local space on a polygon and to
/// evaluate the velocity basis functions, their divergence and the pressure basis functions.
/// Each public method dispatches on @ref FEM_MCC_2D_LocalSpace_Data::fem_type and forwards the
/// call to the concrete element implementation (currently only the Raviart-Thomas triangle,
/// @ref FEM_Triangle_RT_MCC_2D_LocalSpace); an unsupported type raises a std::runtime_error.
class FEM_MCC_2D_LocalSpace final
{
  private:
    Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace rt_triangle_local_space;

  public:
    /// @brief Builds the local space on a given polygon.
    /// @param reference_element_data Reference element data selecting the finite element type and its parameters.
    /// @param polygon Geometric description of the physical polygon on which the space is built.
    /// @return The assembled local space data, including quadrature rules and basis function counts.
    Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data CreateLocalSpace(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry &polygon) const;

    /// @brief Evaluates the (vector-valued) velocity basis functions at the local space internal quadrature points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @return One matrix per spatial component; each matrix holds a basis function per column evaluated at the
    /// quadrature points (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    std::vector<Eigen::MatrixXd> ComputeVelocityBasisFunctionsValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputeVelocityBasisFunctionsValues(reference_element_data.rt_triangle_reference_element_data,
                                                                               local_space.rt_triangle_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    /// @brief Evaluates the (scalar) pressure basis functions at the local space internal quadrature points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @return A matrix holding a basis function per column evaluated at the quadrature points (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    Eigen::MatrixXd ComputePressureBasisFunctionsValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputePressureBasisFunctionsValues(reference_element_data.rt_triangle_reference_element_data,
                                                                               local_space.rt_triangle_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    /// @brief Evaluates the divergence of the velocity basis functions at the local space internal quadrature points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @return A matrix holding the divergence of a basis function per column evaluated at the quadrature points
    /// (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    Eigen::MatrixXd ComputeVelocityBasisFunctionsDivergenceValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                  const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputeVelocityBasisFunctionsDivergenceValues(
                reference_element_data.rt_triangle_reference_element_data,
                local_space.rt_triangle_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    /// @brief Evaluates the (vector-valued) velocity basis functions at a user-provided set of points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @param points Evaluation points, stored one point per column.
    /// @return One matrix per spatial component; each matrix holds a basis function per column evaluated at the given
    /// points (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    std::vector<Eigen::MatrixXd> ComputeVelocityBasisFunctionsValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space,
                                                                     const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputeVelocityBasisFunctionsValues(reference_element_data.rt_triangle_reference_element_data,
                                                                               local_space.rt_triangle_local_space_data,
                                                                               points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    /// @brief Evaluates the (scalar) pressure basis functions at a user-provided set of points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @param points Evaluation points, stored one point per column.
    /// @return A matrix holding a basis function per column evaluated at the given points (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    Eigen::MatrixXd ComputePressureBasisFunctionsValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space,
                                                        const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputePressureBasisFunctionsValues(reference_element_data.rt_triangle_reference_element_data,
                                                                               local_space.rt_triangle_local_space_data,
                                                                               points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    /// @brief Evaluates the divergence of the velocity basis functions at a user-provided set of points.
    /// @param reference_element_data Reference element data for the selected finite element type.
    /// @param local_space Local space data previously built with @ref CreateLocalSpace.
    /// @param points Evaluation points, stored one point per column.
    /// @return A matrix holding the divergence of a basis function per column evaluated at the given points (rows).
    /// @throws std::runtime_error if the finite element type is not supported.
    Eigen::MatrixXd ComputeVelocityBasisFunctionsDivergenceValues(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                  const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data &local_space,
                                                                  const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_2D_Types::RT_Triangle: {

            return rt_triangle_local_space.ComputeVelocityBasisFunctionsDivergenceValues(reference_element_data.rt_triangle_reference_element_data,
                                                                                         local_space.rt_triangle_local_space_data,
                                                                                         points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }
};
} // namespace MCC
} // namespace FEM
} // namespace Polydim

#endif
