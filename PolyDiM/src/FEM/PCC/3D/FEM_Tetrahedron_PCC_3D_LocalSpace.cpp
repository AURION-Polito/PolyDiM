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

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);
    unsigned int dofCounter = 0;
    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 4; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[dofCounter] = d;
            dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[4]);
    for (unsigned int e = 0; e < 6; e++)
    {
        localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

        if (polyhedron.EdgesDirection.at(e))
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e]; d < localSpace.Dof1DsIndex[e + 1]; d++)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e + 1] - 1; d < UINT_MAX && d >= localSpace.Dof1DsIndex[e]; d--)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
    }

    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[6]);
    for (unsigned int f = 0; f < 4; f++)
    {
        localSpace.Dof2DsIndex[f + 1] = localSpace.Dof2DsIndex[f] + reference_element_data.NumDofs2D;

        if (polyhedron.FacesDirection.at(f))
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f]; d--)
            {
                localSpace.DofsMeshOrder[dofCounter] = d;
                dofCounter++;
            }
        }
    }

    localSpace.Dof3DsIndex.fill(localSpace.Dof2DsIndex[4]);
    localSpace.Dof3DsIndex[1] = localSpace.Dof3DsIndex[0] + reference_element_data.NumDofs3D;
    for (unsigned int d = localSpace.Dof3DsIndex[0]; d < localSpace.Dof3DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

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
