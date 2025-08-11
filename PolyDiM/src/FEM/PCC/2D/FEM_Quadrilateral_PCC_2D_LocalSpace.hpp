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

#ifndef __FEM_Quadrilateral_PCC_2D_LocalSpace_HPP
#define __FEM_Quadrilateral_PCC_2D_LocalSpace_HPP

#include "FEM_PCC_2D_LocalSpace_Data.hpp"
#include "FEM_Quadrilateral_PCC_2D_ReferenceElement.hpp"
#include "MapParallelogram.hpp"
#include "MapQuadrilateral.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_Quadrilateral_PCC_2D_LocalSpace final
{
  private:
    Eigen::MatrixXd MapValues(const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

    std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues,
                                                     const Eigen::MatrixXd &referencePoints) const;

    Eigen::MatrixXd MapLaplacianValues(const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                       const std::array<Eigen::MatrixXd, 4> &referenceSecondDerivateValues,
                                       const Eigen::MatrixXd &referencePoints) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space) const;

    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                      const FEM_PCC_2D_Polygon_Geometry &polygon) const;

  public:
    FEM_Quadrilateral_PCC_2D_LocalSpace_Data CreateLocalSpace(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                              const FEM_PCC_2D_Polygon_Geometry &polygon) const;

    Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space) const
    {
        return MapLaplacianValues(local_space,
                                  reference_element_data.ReferenceBasisFunctionSecondDerivativeValues,
                                  reference_element_data.ReferenceSquareQuadrature.Points);
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space,
                                   reference_element_data.ReferenceBasisFunctionDerivativeValues,
                                   reference_element_data.ReferenceSquareQuadrature.Points);
    }

    Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &points) const
    {

        Eigen::MatrixXd referencePoints;
        switch (local_space.quadrilateral_type)
        {
        case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
            Gedim::MapParallelogram mapQuadrilateral;
            referencePoints = mapQuadrilateral.FInv(local_space.MapData, points);
        }
        break;
        case Polydim::FEM::PCC::QuadrilateralType::Generic: {
            Gedim::MapQuadrilateral mapQuadrilateral;
            referencePoints = mapQuadrilateral.FInv(local_space.Vertices, points);
        }
        break;
        default:
            throw std::runtime_error("not valid quadrilateral");
        }

        FEM_Quadrilateral_PCC_2D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
                                                                       const Eigen::MatrixXd &points) const
    {
        Eigen::MatrixXd referencePoints;
        switch (local_space.quadrilateral_type)
        {
        case Polydim::FEM::PCC::QuadrilateralType::Parallelogram: {
            Gedim::MapParallelogram mapQuadrilateral;
            referencePoints = mapQuadrilateral.FInv(local_space.MapData, points);
        }
        break;
        case Polydim::FEM::PCC::QuadrilateralType::Generic: {
            Gedim::MapQuadrilateral mapQuadrilateral;
            referencePoints = mapQuadrilateral.FInv(local_space.Vertices, points);
        }
        break;
        default:
            throw std::runtime_error("not valid quadrilateral");
        }

        FEM_Quadrilateral_PCC_2D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space,
                                   reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data),
                                   referencePoints);
    }

    Eigen::MatrixXd EdgeDOFsCoordinates(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                        const FEM_Quadrilateral_PCC_2D_LocalSpace_Data &local_space,
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

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
    {
        return reference_element_data.BoundaryReferenceElement_Data.ReferenceBasisFunctionValues;
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                      const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::MatrixXd points = Eigen::MatrixXd(3, pointsCurvilinearCoordinates.size());
        points.row(0) = pointsCurvilinearCoordinates;
        FEM_PCC_1D_ReferenceElement reference_element;
        return reference_element.EvaluateBasisFunctions(points, reference_element_data.BoundaryReferenceElement_Data);
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
