#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_Tetrahedron_PCC_3D_LocalSpace_Data FEM_Tetrahedron_PCC_3D_LocalSpace::CreateLocalSpace(
    const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
    const FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = polyhedron.Tolerance1D;
    geometry_utilities_config.Tolerance2D = polyhedron.Tolerance2D;
    geometry_utilities_config.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MapTetrahedron mapTetrehedron(geometry_utilities);
    localSpace.MapData = mapTetrehedron.Compute(polyhedron.Vertices);

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, Gedim::MapTetrahedron::F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTetrahedronQuadrature, localSpace.MapData);

    return localSpace;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_Tetrahedron_PCC_3D_LocalSpace::MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                             const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
{
    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
        3,
        Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int i = 0; i < 3; i++)
    {
        basisFunctionsDerivativeValues[i] = local_space.MapData.QInv(i, i) * referenceDerivateValues[i];
        for (unsigned int j = 0; j < i; j++)
        {
            basisFunctionsDerivativeValues[i] += local_space.MapData.QInv(j, i) * referenceDerivateValues[j];
            basisFunctionsDerivativeValues[j] += local_space.MapData.QInv(i, j) * referenceDerivateValues[i];
        }
    }

    return basisFunctionsDerivativeValues;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Tetrahedron_PCC_3D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    quadrature.Points = Gedim::MapTetrahedron::F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapTetrahedron::DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
