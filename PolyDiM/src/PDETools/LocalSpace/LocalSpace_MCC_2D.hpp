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

#ifndef __LocalSpace_MCC_2D_H
#define __LocalSpace_MCC_2D_H

#include "FEM_MCC_2D_LocalSpace.hpp"
#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_MCC_2D_Creator.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace PDETools
{

/// @brief Unified interface to the 2D mixed conforming (MCC) local spaces.
///
/// This namespace provides a single, method-agnostic API for building and evaluating
/// the coupled velocity(flux)/pressure local spaces of a 2D mixed formulation (e.g. of
/// an elliptic/Darcy problem). It hides the differences between the mixed virtual
/// element method (VEM MCC), in its several orthogonalized/partial variants, and the
/// Raviart–Thomas finite element (FEM RT). The active backend is selected through
/// MethodTypes; the free functions dispatch accordingly, so assemblers can be written
/// once and reused across methods.
namespace LocalSpace_MCC_2D
{
/// @brief Discretization method for the 2D MCC velocity/pressure pair.
enum class MethodTypes
{
    VEM_MCC = 1,                 ///< Mixed VEM (standard basis). \cite secondMixed
    VEM_MCC_Partial = 2,         ///< Mixed VEM with a partially computed velocity space.
    VEM_MCC_Ortho = 3,           ///< Mixed VEM with an orthonormalized internal basis.
    VEM_MCC_EdgeOrtho = 4,       ///< Mixed VEM with an orthonormalized edge basis. \cite Teora2023
    VEM_MCC_Ortho_EdgeOrtho = 5, ///< Mixed VEM with both internal and edge orthonormalization. \cite Teora2023
    FEM_RT_MCC = 6               ///< Raviart–Thomas finite element.
};

/// @brief Method-specific reference-element data for the 2D MCC local space.
///
/// Holds the selected method and polynomial order together with the velocity and
/// pressure reference-element / local-space handles of the active backend. Only the
/// members belonging to @ref Method_Type are populated.
class ReferenceElement_Data final
{
  public:
    Polydim::PDETools::LocalSpace_MCC_2D::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Velocity_ReferenceElement> VEM_ReferenceElement_Velocity;
    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Pressure_ReferenceElement> VEM_ReferenceElement_Pressure;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data VEM_ReferenceElement_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data VEM_ReferenceElement_Data_Pressure;
    Polydim::VEM::MCC::VEM_MCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Velocity_LocalSpace> VEM_LocalSpace_Velocity;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Pressure_LocalSpace> VEM_LocalSpace_Pressure;

    Polydim::FEM::MCC::FEM_MCC_Types FEM_Type;
    Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace FEM_LocalSpace;
    Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement FEM_ReferenceElement;
    Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
};

/// @brief Per-cell local-space data for the 2D MCC local space.
///
/// Stores the polygon geometry and the computed velocity/pressure local-space data for
/// the active backend; only the members of the selected method are populated.
class LocalSpace_Data final
{
  public:
    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data VEM_LocalSpace_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data VEM_LocalSpace_Data_Pressure;

    Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry FEM_Geometry;
    Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace_Data FEM_LocalSpace_Data;
};

/// @brief Per-cell performance metrics for the 2D MCC local space.
class Performance_Data final
{
  public:
    class Cell2D_Performance final
    {
      public:
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data VEM_Analysis;
    };

    Polydim::PDETools::LocalSpace_MCC_2D::Performance_Data::Cell2D_Performance Performance_Data;
};

/// @brief Create the velocity/pressure reference elements for the chosen method and order.
///
/// Factory that instantiates the backend-specific reference elements and populates the
/// corresponding members of the returned @ref ReferenceElement_Data.
///
/// @param method_type  Discretization method to use.
/// @param method_order Polynomial order of the space.
/// @return The initialized reference-element data.
Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_MCC_2D::MethodTypes &method_type,
                                                                                   const unsigned int method_order);

/// @brief Number of reference-element DOFs per mesh-entity type, for velocity and pressure.
///
/// Returns, for each of the two fields, the number of degrees of freedom associated
/// with the four mesh-entity categories (e.g. vertices, edges, internal, and the
/// remaining entity slot), as used to lay out the DOFs.
///
/// @param reference_element_data Reference-element data describing the method.
/// @return A 2\f$\times\f$4 array of DOF counts (row 0: velocity, row 1: pressure).
std::array<std::array<unsigned int, 4>, 2> ReferenceElementNumDOFs(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data);

/// @brief Build the coupled velocity/pressure local space on a given cell.
///
/// Computes the backend-specific velocity and pressure local-space data on cell
/// @p cell2D_index from the precomputed mesh geometric data.
///
/// @param geometric_tolerance_1D Geometric tolerance for 1D (edge) operations.
/// @param geometric_tolerance_2D Geometric tolerance for 2D (polygon) operations.
/// @param mesh_geometric_data    Precomputed geometric data of the mesh.
/// @param cell2D_index           Index of the 2D cell to process.
/// @param reference_element_data Reference-element data describing the method.
/// @return The computed local-space data for the cell.
Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data CreateLocalSpace(
    const double &geometric_tolerance_1D,
    const double &geometric_tolerance_2D,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const unsigned int cell2D_index,
    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data);

/// @brief Evaluate the (vector) velocity/flux basis functions at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM); defaults to \f$\Pi^0_k\f$.
/// @return One (quadrature points \f$\times\f$ velocity DOFs) matrix per velocity component.
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

/// @brief Evaluate the (scalar) pressure basis functions at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return A (quadrature points \f$\times\f$ pressure DOFs) matrix of pressure basis values.
Eigen::MatrixXd PressureBasisFunctionsValues(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

/// @brief Evaluate the velocity/flux divergence at the internal quadrature points.
///
/// Returns \f$\nabla \cdot \boldsymbol{u}\f$ for each velocity basis function, as used
/// in the divergence (velocity–pressure coupling) term of the mixed formulation.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return A (quadrature points \f$\times\f$ velocity DOFs) matrix of divergence values.
Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                       const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

/// @brief Compute the (VEM) stabilization matrix of the velocity local space.
///
/// Returns the stabilization term complementing the consistency part of the mixed VEM
/// velocity bilinear form; empty/zero for the Raviart–Thomas FEM.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         Projection used to define the stabilization; defaults to \f$\Pi^0_k\f$.
/// @return The local stabilization matrix.
Eigen::MatrixXd StabilizationMatrix(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

/// @brief Coordinates (and weights) of the velocity DOFs on a given edge.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_local_index       Local index of the edge within the cell.
/// @return Quadrature-like data holding the edge-DOF coordinates and weights.
Gedim::Quadrature::QuadratureData EdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                      const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                      const unsigned int edge_local_index);

/// @brief Compute the edge velocity DOFs from strongly-imposed boundary values.
///
/// Evaluates the degrees of freedom associated with an edge (the normal-flux moments)
/// from the boundary data sampled at the edge-DOF coordinates, e.g. to strongly impose
/// Neumann/essential conditions of the mixed problem.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_local_index       Local index of the edge within the cell.
/// @param edge_dofs_coordinates  Edge-DOF coordinates/weights from EdgeDofsCoordinates().
/// @param strong_values          Boundary values sampled at the edge-DOF coordinates.
/// @return The vector of edge velocity DOFs.
Eigen::VectorXd EdgeDofs(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                         const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                         const unsigned int edge_local_index,
                         const Gedim::Quadrature::QuadratureData &edge_dofs_coordinates,
                         const Eigen::VectorXd &strong_values);

/// @brief Evaluate the velocity basis functions on a given edge.
///
/// @param edge_local_index       Local index of the edge within the cell.
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_quadrature_points Evaluation points on the edge.
/// @return A matrix of velocity edge basis-function values at the requested points.
Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                                    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                    const Eigen::MatrixXd &edge_quadrature_points);

/// @brief Evaluate the (vector) velocity/flux basis functions at arbitrary points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param points                 Evaluation points (one per column).
/// @param projectionType         VEM projection operator to apply (ignored by FEM); defaults to \f$\Pi^0_k\f$.
/// @return One (points \f$\times\f$ velocity DOFs) matrix per velocity component.
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Eigen::MatrixXd &points,
    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

/// @brief Quadrature rule on a given (physical) edge of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_local_index       Local index of the edge within the cell.
/// @return The edge quadrature points and weights.
Gedim::Quadrature::QuadratureData EdgeQuadrature(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                 const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                 const unsigned int edge_local_index);

/// @brief Quadrature rule on the reference edge.
///
/// @param reference_element_data Reference-element data describing the method.
/// @return The quadrature points and weights on the reference edge, independent of the cell.
Gedim::Quadrature::QuadratureData EdgeReferenceQuadrature(const ReferenceElement_Data &reference_element_data);

/// @brief Internal quadrature rule of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The internal quadrature points and weights used for volume integration.
Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

/// @brief Number of local velocity degrees of freedom of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The velocity local-space dimension (number of local velocity DOFs).
unsigned int VelocitySize(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                          const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

/// @brief Compute the per-cell performance metrics of the local space.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The collected @ref Performance_Data (quadrature counts and VEM analysis).
Polydim::PDETools::LocalSpace_MCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

} // namespace LocalSpace_MCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
