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

#include "FEM_PCC_2D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_PCC_2D_LocalSpace_Data FEM_PCC_2D_LocalSpace::CreateLocalSpace(const FEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const FEM_PCC_2D_Polygon_Geometry &polygon) const
{
    FEM_PCC_2D_LocalSpace_Data localSpace;

    switch (polygon.Vertices.cols())
    {
    case 3: {
        localSpace.fem_type = FEM_PCC_2D_Types::Triangle;
        localSpace.triangle_local_space_data =
            triangle_local_space.CreateLocalSpace(reference_element_data.triangle_reference_element_data, polygon);

        localSpace.InternalQuadrature = localSpace.triangle_local_space_data.InternalQuadrature;
        localSpace.BoundaryQuadrature = localSpace.triangle_local_space_data.BoundaryQuadrature;
        localSpace.NumberOfBasisFunctions = localSpace.triangle_local_space_data.NumberOfBasisFunctions;
    }
    break;
    case 4: {
        localSpace.fem_type = FEM_PCC_2D_Types::Quadrilateral;
        localSpace.quadrilateral_local_space_data =
            quadrilateral_local_space.CreateLocalSpace(reference_element_data.quadrilateral_reference_element_data, polygon);

        localSpace.InternalQuadrature = localSpace.quadrilateral_local_space_data.InternalQuadrature;
        localSpace.BoundaryQuadrature = localSpace.quadrilateral_local_space_data.BoundaryQuadrature;
        localSpace.NumberOfBasisFunctions = localSpace.quadrilateral_local_space_data.NumberOfBasisFunctions;
    }
    break;
    default:
        throw std::runtime_error("not valid element");
    }

    return localSpace;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
