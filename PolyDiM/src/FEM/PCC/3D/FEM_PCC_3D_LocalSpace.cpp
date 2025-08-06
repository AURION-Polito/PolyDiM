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

#include "FEM_PCC_3D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_PCC_3D_LocalSpace_Data FEM_PCC_3D_LocalSpace::CreateLocalSpace(const FEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                   const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    FEM_PCC_3D_LocalSpace_Data localSpace;

    if (polyhedron.Vertices.cols() == 4 && polyhedron.Edges.cols() == 6 && polyhedron.Faces.size() == 4)
    {
        localSpace.fem_type = FEM_PCC_3D_Types::Tetrahedron;
        localSpace.tetrahedron_local_space_data =
            tetrahedron_local_space.CreateLocalSpace(reference_element_data.tetrahedron_reference_element_data, polyhedron);

        localSpace.InternalQuadrature = localSpace.tetrahedron_local_space_data.InternalQuadrature;
        localSpace.BoundaryQuadrature.resize(4);
        for (unsigned int f = 0; f < 4; f++)
            localSpace.BoundaryQuadrature[f] = localSpace.tetrahedron_local_space_data.BoundaryQuadrature[f];
        localSpace.NumberOfBasisFunctions = localSpace.tetrahedron_local_space_data.NumberOfBasisFunctions;
    }
    else if (polyhedron.Vertices.cols() == 8 && polyhedron.Edges.cols() == 12 && polyhedron.Faces.size() == 6)
    {
        localSpace.fem_type = FEM_PCC_3D_Types::Hexahedron;
        localSpace.hexahedron_local_space_data =
            hexahedron_local_space.CreateLocalSpace(reference_element_data.hexahedron_reference_element_data, polyhedron);

        localSpace.InternalQuadrature = localSpace.hexahedron_local_space_data.InternalQuadrature;
        localSpace.BoundaryQuadrature.resize(6);
        for (unsigned int f = 0; f < 6; f++)
            localSpace.BoundaryQuadrature[f] = localSpace.hexahedron_local_space_data.BoundaryQuadrature[f];
        localSpace.NumberOfBasisFunctions = localSpace.hexahedron_local_space_data.NumberOfBasisFunctions;
    }
    else
        throw std::runtime_error("not valid element");

    return localSpace;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
