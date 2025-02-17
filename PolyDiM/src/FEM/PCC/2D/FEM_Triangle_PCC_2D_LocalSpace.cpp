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

    if (polygon.Vertices.cols() != 3)
        throw std::runtime_error("The element must be a triangle");

    FEM_Triangle_PCC_2D_LocalSpace_Data localSpace;

    Gedim::MapTriangle mapTriangle;
    localSpace.MapData = mapTriangle.Compute(polygon.Vertices);
    localSpace.B_lap = localSpace.MapData.BInv * localSpace.MapData.BInv.transpose();

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);
    unsigned int dofCounter = 0;
    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 3; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[dofCounter] = d;
            dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[3]);
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

    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[3]);
    localSpace.Dof2DsIndex[1] = localSpace.Dof2DsIndex[0] + reference_element_data.NumDofs2D;
    for (unsigned int d = localSpace.Dof2DsIndex[0]; d < localSpace.Dof2DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[dofCounter] = d;
        dofCounter++;
    }

    // reorder basis function values with mesh order
    localSpace.Dofs = MapValues(localSpace, Gedim::MapTriangle::F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTriangleQuadrature, localSpace.MapData);
    localSpace.BoundaryQuadrature =
        BoundaryQuadrature(reference_element_data.BoundaryReferenceElement_Data.ReferenceSegmentQuadrature, polygon);

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
std::vector<Gedim::Quadrature::QuadratureData> FEM_Triangle_PCC_2D_LocalSpace::BoundaryQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_Triangle_PCC_2D_Polygon_Geometry &polygon) const
{
    const unsigned int num_edges = polygon.EdgesDirection.size();
    std::vector<Gedim::Quadrature::QuadratureData> edges_quadrature(num_edges);

    for (unsigned int e = 0; e < num_edges; ++e)
    {
        auto &edge_quadrature = edges_quadrature.at(e);
        const bool edge_direction = polygon.EdgesDirection.at(e);
        const double edge_length = polygon.EdgesLength[e];
        const Eigen::Vector3d edge_origin = edge_direction ? polygon.Vertices.col(e) : polygon.Vertices.col((e + 1) % num_edges);

        const Eigen::Vector3d edge_tangent = edge_direction ? +1.0 * polygon.EdgesTangent.col(e)
                                                            : -1.0 * polygon.EdgesTangent.col(e);

        const unsigned int num_quadrature_points = reference_quadrature.Points.cols();
        edge_quadrature.Points.resize(3, num_quadrature_points);
        for (unsigned int q = 0; q < num_quadrature_points; q++)
            edge_quadrature.Points.col(q) = edge_origin + reference_quadrature.Points(0, q) * edge_tangent;

        edge_quadrature.Weights = reference_quadrature.Weights * std::abs(edge_length);
    }

    return edges_quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
