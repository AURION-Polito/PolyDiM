#include "FEM_PCC_1D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_PCC_1D_LocalSpace_Data FEM_PCC_1D_LocalSpace::CreateLocalSpace(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                                   const FEM_PCC_1D_Segment_Geometry &segment) const
{
    FEM_PCC_1D_LocalSpace_Data localSpace;

    localSpace.MapData.Origin = segment.Origin;
    localSpace.MapData.Tangent = segment.Tangent;
    localSpace.MapData.Length = segment.Length;
    localSpace.MapData.SquaredLength = segment.Length * segment.Length;

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceSegmentQuadrature, localSpace.MapData);

    return localSpace;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_PCC_1D_LocalSpace::InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                            const FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    quadrature.Points = F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() * DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
