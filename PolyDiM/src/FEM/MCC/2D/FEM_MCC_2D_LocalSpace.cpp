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

#include "FEM_MCC_2D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace MCC
{
// ***************************************************************************
FEM_MCC_2D_LocalSpace_Data FEM_MCC_2D_LocalSpace::CreateLocalSpace(const FEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const FEM_MCC_2D_Polygon_Geometry &polygon,
                                                                   const Polydim::FEM::MCC::FEM_MCC_Types &fem_main_type) const
{
    FEM_MCC_2D_LocalSpace_Data localSpace;

    switch (polygon.Vertices.cols())
    {
    case 3: {
        switch (fem_main_type)
        {
        case Polydim::FEM::MCC::FEM_MCC_Types::RT: {
            localSpace.fem_type = FEM_MCC_2D_Types::RT_Triangle;
            localSpace.rt_triangle_local_space_data =
                rt_triangle_local_space.CreateLocalSpace(reference_element_data.rt_triangle_reference_element_data, polygon);

            localSpace.InternalQuadrature = localSpace.rt_triangle_local_space_data.InternalQuadrature;
            localSpace.BoundaryQuadrature = localSpace.rt_triangle_local_space_data.BoundaryQuadrature;
            localSpace.NumVelocityBasisFunctions = localSpace.rt_triangle_local_space_data.NumVelocityBasisFunctions;
            localSpace.NumPressureBasisFunctions = localSpace.rt_triangle_local_space_data.NumPressureBasisFunctions;
        }
        break;
        default:
            throw std::runtime_error("not valid mcc fem type");
        }
    }
    break;
    default:
        throw std::runtime_error("not valid element");
    }

    return localSpace;
}
// ***************************************************************************
} // namespace MCC
} // namespace FEM
} // namespace Polydim
