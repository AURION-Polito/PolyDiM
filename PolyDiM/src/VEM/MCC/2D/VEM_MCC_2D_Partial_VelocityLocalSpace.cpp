#include "VEM_MCC_2D_Partial_VelocityLocalSpace.hpp"


#include "LAPACK_utilities.hpp"

using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace MCC
{
//****************************************************************************
VEM_MCC_VelocityLocalSpace_Data VEM_MCC_2D_Partial_VelocityLocalSpace::CreateLocalSpace(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                                                        const VEM_MCC_2D_Polygon_Geometry& polygon) const
{
    VEM_MCC_VelocityLocalSpace_Data localSpace;

    Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature = quadrature.PolygonInternalQuadrature(reference_element_data.Quadrature,
                                                                         polygon.TriangulationVertices);
    localSpace.BoundaryQuadrature = quadrature.PolygonEdgesLobattoQuadrature(reference_element_data.Quadrature,
                                                                             polygon.Vertices,
                                                                             polygon.EdgesLength,
                                                                             polygon.EdgesDirection,
                                                                             polygon.EdgesTangent,
                                                                             polygon.EdgesNormal);


    ComputePolynomialBasisDofs(polygon.Measure,
                               localSpace);

    ComputeStabilizationMatrix(polygon.Measure,
                               localSpace);

    return localSpace;
}
//****************************************************************************

}
}
}
