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

#ifndef __FEM_PCC_3D_LocalSpace_HPP
#define __FEM_PCC_3D_LocalSpace_HPP

#include "FEM_Hexahedron_PCC_3D_LocalSpace.hpp"
#include "FEM_PCC_3D_LocalSpace_Data.hpp"
#include "FEM_PCC_3D_ReferenceElement.hpp"
#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

class FEM_PCC_3D_LocalSpace final
{
private:
    FEM_Tetrahedron_PCC_3D_LocalSpace tetrahedron_local_space;
    FEM_Hexahedron_PCC_3D_LocalSpace hexahedron_local_space;

public:
    FEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_PCC_3D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsValues(reference_element_data.tetrahedron_reference_element_data,
                                                                       local_space.tetrahedron_local_space_data);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsValues(reference_element_data.hexahedron_reference_element_data,
                                                                      local_space.hexahedron_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_PCC_3D_LocalSpace_Data &local_space) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.tetrahedron_reference_element_data,
                                                                                 local_space.tetrahedron_local_space_data);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.hexahedron_reference_element_data,
                                                                                local_space.hexahedron_local_space_data);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_PCC_3D_LocalSpace_Data &local_space,
                                                       const Eigen::MatrixXd &points) const
    {

        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsValues(reference_element_data.tetrahedron_reference_element_data,
                                                                       local_space.tetrahedron_local_space_data,
                                                                       points);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsValues(reference_element_data.hexahedron_reference_element_data,
                                                                      local_space.hexahedron_local_space_data,
                                                                      points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_PCC_3D_LocalSpace_Data &local_space,
                                                                              const Eigen::MatrixXd &points) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.tetrahedron_reference_element_data,
                                                                                 local_space.tetrahedron_local_space_data,
                                                                                 points);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.hexahedron_reference_element_data,
                                                                                local_space.hexahedron_local_space_data,
                                                                                points);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_PCC_3D_LocalSpace_Data &local_space,
                                                             const unsigned int face_index) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsValuesOnFace(reference_element_data.tetrahedron_reference_element_data,
                                                                             local_space.tetrahedron_local_space_data,
                                                                             face_index);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsValuesOnFace(reference_element_data.hexahedron_reference_element_data,
                                                                            local_space.hexahedron_local_space_data,
                                                                            face_index);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_PCC_3D_LocalSpace_Data &local_space,
                                                             const unsigned int face_index,
                                                             const Eigen::MatrixXd &points2D) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.ComputeBasisFunctionsValuesOnFace(reference_element_data.tetrahedron_reference_element_data,
                                                                             local_space.tetrahedron_local_space_data,
                                                                             face_index,
                                                                             points2D);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsValuesOnFace(reference_element_data.hexahedron_reference_element_data,
                                                                            local_space.hexahedron_local_space_data,
                                                                            face_index,
                                                                            points2D);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd FaceDOFsCoordinates(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                               const FEM_PCC_3D_LocalSpace_Data &local_space,
                                               const unsigned int face_local_index) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.FaceDOFsCoordinates(reference_element_data.tetrahedron_reference_element_data,
                                                               local_space.tetrahedron_local_space_data,
                                                               face_local_index);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.FaceDOFsCoordinates(reference_element_data.hexahedron_reference_element_data,
                                                              local_space.hexahedron_local_space_data,
                                                              face_local_index);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd EdgeDOFsCoordinates(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                               const FEM_PCC_3D_LocalSpace_Data &local_space,
                                               const unsigned int edge_local_index) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {

            return tetrahedron_local_space.EdgeDOFsCoordinates(reference_element_data.tetrahedron_reference_element_data,
                                                               local_space.tetrahedron_local_space_data,
                                                               edge_local_index);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.EdgeDOFsCoordinates(reference_element_data.hexahedron_reference_element_data,
                                                              local_space.hexahedron_local_space_data,
                                                              edge_local_index);
        }
        default:
            throw std::runtime_error("not valid fem type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_PCC_3D_LocalSpace_Data &local_space,
                                                             const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        switch (local_space.fem_type)
        {
        case FEM_PCC_3D_Types::Tetrahedron: {
            return tetrahedron_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.tetrahedron_reference_element_data,
                                                                             pointsCurvilinearCoordinates);
        }
        case FEM_PCC_3D_Types::Hexahedron: {
            return hexahedron_local_space.ComputeBasisFunctionsValuesOnEdge(reference_element_data.hexahedron_reference_element_data,
                                                                            pointsCurvilinearCoordinates);
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
