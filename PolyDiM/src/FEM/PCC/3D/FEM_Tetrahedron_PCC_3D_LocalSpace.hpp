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

#ifndef __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP
#define __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP

#include "FEM_PCC_3D_LocalSpace_Data.hpp"
#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"
#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "MapTetrahedron.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_Tetrahedron_PCC_3D_LocalSpace final
{
  private:
    Eigen::MatrixXd MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

    std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const;
    std::array<Gedim::Quadrature::QuadratureData, 4> BoundaryQuadrature(const std::array<FEM_Triangle_PCC_2D_LocalSpace_Data, 4> &faces_local_space_data,
                                                                        const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

  public:
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data CreateLocalSpace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                            const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = Gedim::MapTetrahedron::FInv(local_space.MapData, points);

        FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                       const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = Gedim::MapTetrahedron::FInv(local_space.MapData, points);

        FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }

    inline Eigen::MatrixXd EdgeDOFsCoordinates(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                        const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
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

    inline Eigen::MatrixXd FaceDOFsCoordinates(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                        const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                        const unsigned int face_local_index) const
    {
        const auto &dof_coordinates = local_space.Dofs;

        const unsigned int cell2DStartingLocalIdex = local_space.Dof2DsIndex.at(face_local_index);
        const unsigned int num_face_dofs = reference_element_data.NumDofs2D;

        return (num_face_dofs == 0) ? Eigen::MatrixXd(0, 0) : dof_coordinates.block(0, cell2DStartingLocalIdex, 3, num_face_dofs);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                      const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                      const unsigned int face_index) const
    {
        FEM_Triangle_PCC_2D_LocalSpace face_local_space;

        return face_local_space.ComputeBasisFunctionsValues(reference_element_data.BoundaryReferenceElement_Data,
                                                            local_space.Boundary_LocalSpace_Data[face_index]);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                             const unsigned int face_index,
                                                             const Eigen::MatrixXd &points2D) const
    {
        FEM_Triangle_PCC_2D_LocalSpace face_local_space;

        return face_local_space.ComputeBasisFunctionsValues(reference_element_data.BoundaryReferenceElement_Data,
                                                            local_space.Boundary_LocalSpace_Data[face_index],
                                                            points2D);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                      const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::MatrixXd points = Eigen::MatrixXd(3, pointsCurvilinearCoordinates.size());
        points.row(0) = pointsCurvilinearCoordinates;
        FEM_PCC_1D_ReferenceElement reference_element;
        return reference_element.EvaluateBasisFunctions(points, reference_element_data.EdgeReferenceElement_Data);
    }


};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
