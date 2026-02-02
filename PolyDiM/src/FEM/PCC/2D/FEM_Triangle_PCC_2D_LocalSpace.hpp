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

#ifndef __FEM_Triangle_PCC_2D_LocalSpace_HPP
#define __FEM_Triangle_PCC_2D_LocalSpace_HPP

#include "FEM_PCC_2D_LocalSpace_Data.hpp"
#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "MapTriangle.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_Triangle_PCC_2D_LocalSpace final
{
  private:
    inline Eigen::MatrixXd MapValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                     const Eigen::MatrixXd &referenceValues) const
    {
        return referenceValues * local_space.dofs_permutation;
    }

    std::vector<Eigen::MatrixXd> MapDerivativeValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Eigen::MatrixXd MapLaplacianValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                       const std::array<Eigen::MatrixXd, 4> &referenceSecondDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Gedim::MapTriangle::MapTriangleData &mapData) const;

    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                      const Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry &polygon) const;

  public:
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data CreateLocalSpace(
        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry &polygon) const;

    Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    Eigen::MatrixXd EdgeDOFsCoordinates(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                        const unsigned int edge_local_index) const
    {
        const auto &dof_coordinates = local_space.Dofs;

        const unsigned int cell1DStartingLocalIdex = local_space.Dof1DsIndex.at(edge_local_index);
        const unsigned int num_edge_dofs = reference_element_data.NumDofs1D;

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const Eigen::MatrixXd edge_dofs_coordinates = dof_coordinates.block(0, cell1DStartingLocalIdex, 3, num_edge_dofs);

        return edge_dofs_coordinates;
    }

    Eigen::MatrixXd InternalDOFsCoordinates(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                            const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
    {
        const auto &dof_coordinates = local_space.Dofs;

        const unsigned int starting_index = local_space.Dof2DsIndex.at(0);
        const unsigned int num_internal_dofs = reference_element_data.NumDofs2D;

        if (num_internal_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const Eigen::MatrixXd face_dofs_coordinates = dof_coordinates.block(0, starting_index, 3, num_internal_dofs);

        return face_dofs_coordinates;
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
    {
        if (local_space.Order > 2)
            throw std::runtime_error("Unsupported order");

        return MapLaplacianValues(local_space, reference_element_data.ReferenceBasisFunctionSecondDerivativeValues);
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &points) const
    {
        Gedim::MapTriangle mapTriangle;
        const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

        Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
        const Eigen::MatrixXd &points) const
    {
        Gedim::MapTriangle mapTriangle;
        const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

        Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data) const
    {
        return reference_element_data.BoundaryReferenceElement_Data.ReferenceBasisFunctionValues;
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                      const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::MatrixXd points = Eigen::MatrixXd(3, pointsCurvilinearCoordinates.size());
        points.row(0) = pointsCurvilinearCoordinates;
        Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement reference_element;
        return reference_element.EvaluateBasisFunctions(points, reference_element_data.BoundaryReferenceElement_Data);
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
