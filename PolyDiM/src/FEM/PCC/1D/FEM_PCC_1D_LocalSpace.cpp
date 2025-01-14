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

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);
    unsigned int dofCounter = 0;
    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 2; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[dofCounter] = d;
            dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[2]);
    localSpace.Dof1DsIndex[1] = localSpace.Dof1DsIndex[0] + reference_element_data.NumDofs1D;
    for (unsigned int d = localSpace.Dof1DsIndex[0]; d < localSpace.Dof1DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceSegmentQuadrature, localSpace.MapData);

    return localSpace;
}
// ***************************************************************************
MatrixXd FEM_PCC_1D_LocalSpace::MapValues(const FEM_PCC_1D_LocalSpace_Data& local_space, const Eigen::MatrixXd& referenceValues) const
{
  Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

  for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
      basisFunctionValuesOrdered.col(local_space.DofsMeshOrder.at(d)) << referenceValues.col(d);

  return basisFunctionValuesOrdered;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_PCC_1D_LocalSpace::MapDerivativeValues(const FEM_PCC_1D_LocalSpace_Data& local_space, const std::vector<Eigen::MatrixXd>& referenceDerivateValues) const
{
  std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(1);

  basisFunctionsDerivativeValues[0] = referenceDerivateValues[0] / local_space.MapData.Length;

  std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
      1,
      Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

  for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
  {
      const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
      basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
  }

  return basisFunctionsDerivativeValuesOrdered;
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
