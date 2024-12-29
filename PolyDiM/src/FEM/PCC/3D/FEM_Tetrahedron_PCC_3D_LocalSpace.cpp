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

    Gedim::MapTriangle mapTriangle;
    localSpace.MapData = mapTriangle.Compute(polyhedron.Vertices);

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);
    unsigned int dofCounter = 0;
    localSpace.Dof0DsIndex.resize(4, 0);
    for (unsigned int v = 0; v < 3; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[dofCounter] = d;
            dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.resize(4, localSpace.Dof0DsIndex[3]);

    localSpace.Dof2DsIndex.resize(2, localSpace.Dof1DsIndex[3]);
    localSpace.Dof2DsIndex[1] = localSpace.Dof2DsIndex[0] + reference_element_data.NumDofs2D;
    for (unsigned int d = localSpace.Dof2DsIndex[0]; d < localSpace.Dof2DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

    localSpace.Dof3DsIndex.resize(0, 0);

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, mapTriangle.F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTriangleQuadrature, localSpace.MapData);

    return localSpace;
}
// ***************************************************************************
MatrixXd FEM_Tetrahedron_PCC_3D_LocalSpace::MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                      const Eigen::MatrixXd &referenceValues) const
{
    Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
    {
      const unsigned int dofOrder = d;
        basisFunctionValuesOrdered.col(dofOrder) << referenceValues.col(d);
    }

    return basisFunctionValuesOrdered;
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

    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
        3,
        Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
    {
        const unsigned int dofOrder = d;
        basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
        basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
    }

    return basisFunctionsDerivativeValuesOrdered;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Tetrahedron_PCC_3D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const Gedim::MapTriangle::MapTriangleData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    Gedim::MapTriangle mapTriangle;
    quadrature.Points = mapTriangle.F(mapData, reference_quadrature.Points);
    quadrature.Weights =
        reference_quadrature.Weights.array() * mapTriangle.DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
