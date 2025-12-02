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

class FEM_RT_MCC_2D_LocalSpace final
{
  private:
    Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace rt_triangle_local_space;

  public:
    Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data CreateLocalSpace(const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry &polygon,
                                                                   const Polydim::FEM::MCC::FEM_MCC_Types &fem_main_type) const;

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
