#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_Triangle_PCC_2D_LocalSpace_Data FEM_Triangle_PCC_2D_LocalSpace::CreateLocalSpace(
    const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
    const FEM_Triangle_PCC_2D_Polygon_Geometry &polygon) const
{
    FEM_Triangle_PCC_2D_LocalSpace_Data localSpace;

    Gedim::MapTriangle mapTriangle;
    localSpace.MapData = mapTriangle.Compute(polygon.Vertices);
    localSpace.B_lap = localSpace.MapData.BInv * localSpace.MapData.BInv.transpose();

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
    for (unsigned int e = 0; e < 3; e++)
    {
        localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

        if (polygon.EdgesDirection.at(e))
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

    localSpace.Dof2DsIndex.resize(2, localSpace.Dof1DsIndex[3]);
    localSpace.Dof2DsIndex[1] = localSpace.Dof2DsIndex[0] + reference_element_data.NumDofs2D;
    for (unsigned int d = localSpace.Dof2DsIndex[0]; d < localSpace.Dof2DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

    localSpace.Dof3DsIndex.resize(0, 0);

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, Gedim::MapTriangle::F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTriangleQuadrature, localSpace.MapData);

    return localSpace;
}
// ***************************************************************************
MatrixXd FEM_Triangle_PCC_2D_LocalSpace::MapValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                   const Eigen::MatrixXd &referenceValues) const
{
    Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
        basisFunctionValuesOrdered.col(local_space.DofsMeshOrder.at(d)) << referenceValues.col(d);

    return basisFunctionValuesOrdered;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_Triangle_PCC_2D_LocalSpace::MapDerivativeValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                                          const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
{
    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
        2,
        Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int i = 0; i < 2; i++)
    {
        basisFunctionsDerivativeValues[i] = local_space.MapData.BInv(i, i) * referenceDerivateValues[i];
        for (unsigned int j = 0; j < i; j++)
        {
            basisFunctionsDerivativeValues[i] += local_space.MapData.BInv(j, i) * referenceDerivateValues[j];
            basisFunctionsDerivativeValues[j] += local_space.MapData.BInv(i, j) * referenceDerivateValues[i];
        }
    }

    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
        2,
        Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
    {
        const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
        basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
        basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
    }

    return basisFunctionsDerivativeValuesOrdered;
}
// ***************************************************************************
Eigen::MatrixXd FEM_Triangle_PCC_2D_LocalSpace::MapLaplacianValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                                   const std::array<Eigen::MatrixXd, 4> &referenceSecondDerivateValues) const
{
    Eigen::MatrixXd laplacian_values =
        Eigen::MatrixXd::Zero(referenceSecondDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions);

    laplacian_values = local_space.B_lap(0, 0) * referenceSecondDerivateValues.at(0) +
                       local_space.B_lap(0, 1) * referenceSecondDerivateValues.at(1) +
                       local_space.B_lap(1, 0) * referenceSecondDerivateValues.at(2) +
                       local_space.B_lap(1, 1) * referenceSecondDerivateValues.at(3);

    return laplacian_values;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Triangle_PCC_2D_LocalSpace::InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                                     const Gedim::MapTriangle::MapTriangleData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    Gedim::MapTriangle mapTriangle;
    quadrature.Points = Gedim::MapTriangle::F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapTriangle::DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
