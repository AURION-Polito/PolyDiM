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

#include "ZFEM_PCC_2D_LocalSpace.hpp"
#include "GeometryUtilities.hpp"
#include "VEM_Quadrature_2D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{
//****************************************************************************
ZFEM_PCC_2D_LocalSpace_Data ZFEM_PCC_2D_LocalSpace::CreateLocalSpace(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const ZFEM_PCC_2D_Polygon_Geometry &polygon) const
{
    ZFEM_PCC_2D_LocalSpace_Data localSpace;
    localSpace.Dimension = 2;
    localSpace.Order = reference_element_data.Order;
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    localSpace.KernelIncenter = polygon.ChebishevCenter;
    localSpace.KernelInRadius = polygon.InRadius;

    const unsigned int num_vertices = polygon.Vertices.cols();

    localSpace.fem_geometry.resize(num_vertices);
    localSpace.fem_local_space_data.resize(num_vertices);

    localSpace.NumVertexBasisFunctions = num_vertices * reference_element_data.NumDofs0D;
    localSpace.NumEdgeBasisFunctions = num_vertices * reference_element_data.NumDofs1D;
    localSpace.NumBoundaryBasisFunctions = localSpace.NumVertexBasisFunctions + localSpace.NumEdgeBasisFunctions;
    localSpace.NumInternalBasisFunctions = reference_element_data.NumDofs2D;
    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions + localSpace.NumInternalBasisFunctions;

    localSpace.NumTotalBasisFunctions = localSpace.NumBoundaryBasisFunctions + 1 +
                                        num_vertices * (reference_element_data.NumDofs1D + reference_element_data.NumDofs2D);
    localSpace.NumVirtualBasisFunctions = localSpace.NumTotalBasisFunctions - localSpace.NumBasisFunctions;

    localSpace.ReferenceEdgeDOFsInternalPoints =
        reference_element_data.fem_reference_element_data.BoundaryReferenceElement_Data.DofPositions.row(0).segment(
            2,
            reference_element_data.Order - 1);

    localSpace.EdgeBasisCoefficients =
        ZFEM_utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, localSpace.ReferenceEdgeDOFsInternalPoints);

    VEM::Quadrature::VEM_Quadrature_2D quadrature;
    localSpace.InternalQuadrature =
        quadrature.PolygonInternalQuadrature(reference_element_data.fem_reference_element_data.ReferenceTriangleQuadrature,
                                             polygon.Traingulations);

    std::vector<Eigen::MatrixXd> triangles_internal_quadrature_points(num_vertices);
    for (unsigned int v = 0; v < num_vertices; v++)
    {
        const Eigen::MatrixXd triangle_vertices = polygon.Traingulations[v];
        const auto triangle_edge_lengths = geometryUtilities.PolygonEdgeLengths(triangle_vertices);
        const auto triangle_edge_tangents = geometryUtilities.PolygonEdgeTangents(triangle_vertices);
        const std::vector<bool> triangle_edge_directions = {true, polygon.EdgesDirection[v], false};

        localSpace.fem_geometry[v] = {polygon.Tolerance1D, polygon.Tolerance2D, triangle_vertices, triangle_edge_directions, triangle_edge_tangents, triangle_edge_lengths};

        localSpace.fem_local_space_data[v] =
            fem_local_space.CreateLocalSpace(reference_element_data.fem_reference_element_data, localSpace.fem_geometry[v]);

        triangles_internal_quadrature_points[v] = localSpace.fem_local_space_data[v].InternalQuadrature.Points;
    }

    localSpace.local_to_total = ZFEM_utilities.CreateMaps(localSpace.Order,
                                                          num_vertices,
                                                          reference_element_data.fem_reference_element_data.NumDofs1D,
                                                          reference_element_data.fem_reference_element_data.NumDofs2D,
                                                          localSpace.NumBasisFunctions);

    PolygonFineNodes(num_vertices,
                     localSpace.NumBasisFunctions,
                     localSpace.NumTotalBasisFunctions,
                     localSpace.NumVirtualBasisFunctions,
                     localSpace.fem_local_space_data,
                     localSpace.local_to_total,
                     localSpace.DOFsCoordinates,
                     localSpace.VirtualNodes,
                     localSpace.FinerNodes);

    localSpace.Nk = reference_element_data.monomials_data.NumMonomials;

    localSpace.VanderInternal = monomials.Vander(reference_element_data.monomials_data,
                                                 localSpace.InternalQuadrature.Points,
                                                 polygon.ChebishevCenter,
                                                 polygon.Diameter);

    localSpace.VanderInternalDerivatives =
        monomials.VanderDerivatives(reference_element_data.monomials_data, localSpace.VanderInternal, polygon.Diameter);

    ComputePolynomialsDofs(reference_element_data, polygon.ChebishevCenter, polygon.Diameter, localSpace);

    localSpace.fem_basis_functions_values = ZFEM_utilities.ComputeFEMBasisFunctionsValues(reference_element_data,
                                                                                          localSpace.NumBasisFunctions,
                                                                                          localSpace.NumVirtualBasisFunctions,
                                                                                          localSpace.fem_local_space_data,
                                                                                          localSpace.local_to_total,
                                                                                          triangles_internal_quadrature_points);

    localSpace.fem_basis_functions_derivative_values =
        ZFEM_utilities.ComputeFEMBasisFunctionsDerivativeValues(localSpace.Dimension,
                                                                reference_element_data,
                                                                localSpace.NumBasisFunctions,
                                                                localSpace.NumVirtualBasisFunctions,
                                                                localSpace.fem_local_space_data,
                                                                localSpace.local_to_total,
                                                                triangles_internal_quadrature_points);

    ZFEM_utilities.ComputeMinimizerSumOfSquaredWeightsMonomials(localSpace.VirtualNodes,
                                                                localSpace.NumBasisFunctions,
                                                                localSpace.Dmatrix,
                                                                localSpace.VanderVirtuals,
                                                                localSpace.VirtualWeights);

    localSpace.ZFEM_basis_functions_values = ZFEM_utilities.ComputeBasisFunctionsValues(localSpace.NumBasisFunctions,
                                                                                        localSpace.NumVirtualBasisFunctions,
                                                                                        localSpace.fem_basis_functions_values,
                                                                                        localSpace.VirtualWeights);

    localSpace.ZFEM_basis_functions_derivative_values =
        ZFEM_utilities.ComputeBasisFunctionsDerivativeValues(localSpace.Dimension,
                                                             localSpace.NumBasisFunctions,
                                                             localSpace.NumVirtualBasisFunctions,
                                                             localSpace.fem_basis_functions_derivative_values,
                                                             localSpace.VirtualWeights);

    return localSpace;
}
//****************************************************************************
void ZFEM_PCC_2D_LocalSpace::ComputePolynomialsDofs(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const Eigen::Vector3d &internal_points,
                                                    const double &diameter,
                                                    ZFEM_PCC_2D_LocalSpace_Data &localSpace) const
{

    localSpace.Dmatrix = MatrixXd::Zero(localSpace.NumBasisFunctions, localSpace.Nk);

    localSpace.VanderDOFs =
        monomials.Vander(reference_element_data.monomials_data, localSpace.DOFsCoordinates, internal_points, diameter);

    // boundary degrees of freedom of monomials (values at points on
    // the boundary)
    localSpace.Dmatrix = localSpace.VanderDOFs;

    localSpace.VanderVirtuals =
        monomials.Vander(reference_element_data.monomials_data, localSpace.VirtualNodes, internal_points, diameter);
}
//****************************************************************************
void ZFEM_PCC_2D_LocalSpace::PolygonFineNodes(const unsigned int num_vertices,
                                              const unsigned int &NumBasisFunctions,
                                              const unsigned int &NumTotalBasisFunctions,
                                              const unsigned int &NumVirtualBasisFunctions,
                                              const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
                                              const Eigen::MatrixXi &local_to_total,
                                              Eigen::MatrixXd &CoarseNodes,
                                              Eigen::MatrixXd &VirtualNodes,
                                              Eigen::MatrixXd &FinerNodes) const
{

    FinerNodes = Eigen::MatrixXd::Zero(3, NumTotalBasisFunctions);

    for (unsigned int t = 0; t < num_vertices; t++)
        FinerNodes(all, local_to_total.row(t)) = fem_local_space_data[t].Dofs;

    CoarseNodes = FinerNodes.leftCols(NumBasisFunctions);
    VirtualNodes = FinerNodes.rightCols(NumVirtualBasisFunctions);
}
//****************************************************************************
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim
