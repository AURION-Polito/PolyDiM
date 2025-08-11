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

#ifndef __FEM_Hexahedron_PCC_3D_LocalSpace_HPP
#define __FEM_Hexahedron_PCC_3D_LocalSpace_HPP

#include "FEM_Hexahedron_PCC_3D_ReferenceElement.hpp"
#include "FEM_PCC_3D_LocalSpace_Data.hpp"
#include "FEM_Quadrilateral_PCC_2D_LocalSpace.hpp"
#include "MapHexahedron.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_Hexahedron_PCC_3D_LocalSpace final
{
  private:
    Eigen::MatrixXd MapValues(const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

    std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space) const;

    std::array<Gedim::Quadrature::QuadratureData, 6> BoundaryQuadrature(const std::array<FEM_Quadrilateral_PCC_2D_LocalSpace_Data, 6> &faces_local_space_data,
                                                                        const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

  public:
    FEM_Hexahedron_PCC_3D_LocalSpace_Data CreateLocalSpace(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                           const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &points) const
    {
        Eigen::MatrixXd referencePoints;
        switch (local_space.hexahedron_type)
        {
        case Polydim::FEM::PCC::HexahedronType::Parallelepiped: {
            referencePoints = Gedim::MapParallelepiped::FInv(local_space.MapDataParallelepiped, points);
        }
        break;
        default:
            throw std::runtime_error("not valid hexahedron type");
        }

        FEM_Hexahedron_PCC_3D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                       const Eigen::MatrixXd &points) const
    {
        Eigen::MatrixXd referencePoints;
        switch (local_space.hexahedron_type)
        {
        case Polydim::FEM::PCC::HexahedronType::Parallelepiped: {
            referencePoints = Gedim::MapParallelepiped::FInv(local_space.MapDataParallelepiped, points);
        }
        break;
        default:
            throw std::runtime_error("not valid hexahedron type");
        }

        FEM_Hexahedron_PCC_3D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }

    Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                      const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                      const unsigned int face_index) const
    {
        FEM_Quadrilateral_PCC_2D_LocalSpace face_local_space;

        return face_local_space.ComputeBasisFunctionsValues(reference_element_data.BoundaryReferenceElement_Data,
                                                            local_space.Boundary_LocalSpace_Data[face_index]);
    }

    Eigen::MatrixXd EdgeDOFsCoordinates(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                        const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
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

    Eigen::MatrixXd FaceDOFsCoordinates(const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                        const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                        const unsigned int face_local_index) const
    {
        const auto &dof_coordinates = local_space.Dofs;

        const unsigned int cell2DStartingLocalIdex = local_space.Dof2DsIndex.at(face_local_index);
        const unsigned int num_face_dofs = reference_element_data.NumDofs2D;

        return (num_face_dofs == 0) ? Eigen::MatrixXd(0, 0) : dof_coordinates.block(0, cell2DStartingLocalIdex, 3, num_face_dofs);
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
