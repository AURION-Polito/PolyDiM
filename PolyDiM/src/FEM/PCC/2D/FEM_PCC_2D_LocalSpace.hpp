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

#ifndef __FEM_PCC_2D_LocalSpace_HPP
#define __FEM_PCC_2D_LocalSpace_HPP

#include "FEM_PCC_2D_LocalSpace_Data.hpp"
#include "FEM_PCC_2D_ReferenceElement.hpp"
#include "FEM_Quadrilateral_PCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_PCC_2D_LocalSpace final
{
  private:
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace triangle_local_space;
    Polydim::FEM::PCC::FEM_Quadrilateral_PCC_2D_LocalSpace quadrilateral_local_space;

  public:
    Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry &polygon) const;

    Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsValues(reference_element_data.triangle_reference_element_data,
                                                                    local_space.triangle_local_space_data);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsValues(reference_element_data.quadrilateral_reference_element_data,
                                                                         local_space.quadrilateral_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsLaplacianValues(reference_element_data.triangle_reference_element_data,
                                                                             local_space.triangle_local_space_data);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsLaplacianValues(reference_element_data.quadrilateral_reference_element_data,
                                                                                  local_space.quadrilateral_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                       const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.triangle_reference_element_data,
                                                                              local_space.triangle_local_space_data);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.quadrilateral_reference_element_data,
                                                                                   local_space.quadrilateral_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsValues(reference_element_data.triangle_reference_element_data,
                                                                    local_space.triangle_local_space_data,
                                                                    points);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsValues(reference_element_data.quadrilateral_reference_element_data,
                                                                         local_space.quadrilateral_local_space_data,
                                                                         points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                       const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space,
                                                                       const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.triangle_reference_element_data,
                                                                              local_space.triangle_local_space_data,
                                                                              points);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.quadrilateral_reference_element_data,
                                                                                   local_space.quadrilateral_local_space_data,
                                                                                   points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                      const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.triangle_reference_element_data);
        }
        case FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.quadrilateral_reference_element_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                      const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space,
                                                      const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {

            return triangle_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.triangle_reference_element_data,
                                                                          pointsCurvilinearCoordinates);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.quadrilateral_reference_element_data,
                                                                               pointsCurvilinearCoordinates);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    Eigen::MatrixXd EdgeDOFsCoordinates(const Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                        const Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data &local_space,
                                        const unsigned int edge_local_index) const
    {
        switch (local_space.fem_type)
        {
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Triangle: {
            return triangle_local_space.EdgeDOFsCoordinates(reference_element_data.triangle_reference_element_data,
                                                            local_space.triangle_local_space_data,
                                                            edge_local_index);
        }
        case Polydim::FEM::PCC::FEM_PCC_2D_Types::Quadrilateral: {
            return quadrilateral_local_space.EdgeDOFsCoordinates(reference_element_data.quadrilateral_reference_element_data,
                                                                 local_space.quadrilateral_local_space_data,
                                                                 edge_local_index);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
