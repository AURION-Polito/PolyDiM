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

#ifndef __LocalSpace_DF_PCC_2D_H
#define __LocalSpace_DF_PCC_2D_H

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "FEM_PCC_2D_LocalSpace.hpp"
#include "IArray.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_2D_Creator.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_PerformanceAnalysis.hpp"
#include <memory>

namespace Polydim
{
namespace PDETools
{

/// @brief Unified interface to the 2D divergence-free primal conforming (DF_PCC) local spaces.
///
/// This namespace provides a single, method-agnostic API for building and evaluating
/// the coupled velocity/pressure local spaces of a 2D incompressible (divergence-free)
/// vector problem — e.g. Stokes, Navier–Stokes or Brinkman flow. It hides the
/// differences between a Taylor–Hood finite element pair and the divergence-free
/// virtual element method (VEM DF), in its full and reduced variants. The active
/// backend is selected through MethodTypes; the free functions dispatch accordingly,
/// so assemblers can be written once and reused across methods.
namespace LocalSpace_DF_PCC_2D
{
/// @brief Discretization method for the 2D DF_PCC velocity/pressure pair.
enum struct MethodTypes
{
    TAYLOR_HOOD = 0,       ///< Taylor–Hood finite element velocity/pressure pair.
    VEM_DF_PCC_FULL = 1,   ///< Divergence-free VEM, full velocity space. \cite DaVeigaLovadina2017
    VEM_DF_PCC_REDUCED = 2 ///< Divergence-free VEM, reduced velocity space. \cite DaVeigaLovadina2017
};

/// @brief Method-specific reference-element data for the 2D DF_PCC local space.
///
/// Holds the selected method, polynomial order and spatial dimension together with the
/// velocity and pressure reference-element / local-space handles of the active backend.
/// Only the members belonging to @ref Method_Type are populated.
class ReferenceElement_Data final
{
  public:
    Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes Method_Type;
    unsigned int Order;
    unsigned int Dimension = 2;

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_ReferenceElement> VEM_Velocity_ReferenceElement;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data VEM_Velocity_ReferenceElement_Data;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_ReferenceElement> VEM_Pressure_ReferenceElement;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data VEM_Pressure_ReferenceElement_Data;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace> VEM_Velocity_LocalSpace;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace> VEM_Pressure_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement> FEM_ReferenceElement;
    Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace> FEM_LocalSpace;
};

/// @brief Per-cell local-space data for the 2D DF_PCC local space.
///
/// Stores the polygon geometry and the computed velocity/pressure local-space data for
/// the active backend; only the members of the selected method are populated.
class LocalSpace_Data final
{
  public:
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data VEM_Velocity_LocalSpace_Data;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_LocalSpace_Data VEM_Pressure_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data FEM_LocalSpace_Data;
};

/// @brief Per-cell performance metrics for the 2D DF_PCC local space.
class Performance_Data final
{
  public:
    class Cell2D_Performance final
    {
      public:
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis_Data vem_analysis_data;
    };

    Polydim::PDETools::LocalSpace_DF_PCC_2D::Performance_Data::Cell2D_Performance performance_data;
};

/// @brief Create the velocity/pressure reference elements for the chosen method and order.
///
/// Factory that instantiates the backend-specific reference elements and populates the
/// corresponding members of the returned @ref ReferenceElement_Data.
///
/// @param method_type  Discretization method to use.
/// @param method_order Polynomial order of the velocity space.
/// @return The initialized reference-element data.
Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes &method_type,
                                                                                      const unsigned int method_order);

/// @brief Build the DOF layout over the mesh for the velocity and pressure fields.
///
/// Determines how the DOFs of the coupled spaces are distributed over mesh entities and
/// marks boundary DOFs from @p boundary_info. The returned vector holds one
/// MeshDOFsInfo per field (velocity and pressure).
///
/// @param reference_element_data Reference-element data describing the method.
/// @param mesh                   Mesh data access object.
/// @param boundary_info          Per-field boundary maps (index 0: velocity, index 1: pressure).
/// @return One mesh DOFs layout per field.
std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> SetMeshDOFsInfo(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::array<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryMap, 2> &boundary_info);

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
LocalSpace_Data CreateLocalSpace(const double &geometric_tolerance_1D,
                                 const double &geometric_tolerance_2D,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data);

/// @brief Evaluate the (vector) velocity basis functions at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM); defaults to \f$\Pi^0_k\f$.
/// @return One (quadrature points \f$\times\f$ velocity DOFs) matrix per velocity component.
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

/// @brief Evaluate the (vector) velocity basis functions at arbitrary points.
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
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

/// @brief Evaluate the (scalar) pressure basis functions at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return A (quadrature points \f$\times\f$ pressure DOFs) matrix of pressure basis values.
Eigen::MatrixXd PressureBasisFunctionsValues(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Evaluate the (scalar) pressure basis functions at arbitrary points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param points                 Evaluation points (one per column).
/// @return A (points \f$\times\f$ pressure DOFs) matrix of pressure basis values.
Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                             const LocalSpace_Data &local_space_data,
                                             const Eigen::MatrixXd &points);

/// @brief Evaluate the velocity basis-function derivatives at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM); defaults to the energy projection
/// \f$\Pi^\nabla\f$.
/// @return The velocity gradient components (one matrix per derivative/component pair).
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

/// @brief Evaluate the velocity basis-function derivatives at arbitrary points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param points                 Evaluation points (one per column).
/// @param projectionType         VEM projection operator to apply (ignored by FEM); defaults to the energy projection
/// \f$\Pi^\nabla\f$.
/// @return The velocity gradient components (one matrix per derivative/component pair).
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Eigen::MatrixXd &points,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

/// @brief Evaluate the velocity divergence at the internal quadrature points.
///
/// Returns \f$\nabla \cdot \boldsymbol{u}\f$ for each velocity basis function; for the
/// divergence-free VEM this quantity is exactly represented.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return A (quadrature points \f$\times\f$ velocity DOFs) matrix of divergence values.
Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                       const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Internal quadrature rule of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The internal quadrature points and weights used for volume integration.
Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Evaluate the velocity basis functions on a given edge.
///
/// @param edge_local_index             Local index of the edge within the cell.
/// @param reference_element_data       Reference-element data describing the method.
/// @param local_space_data             Local-space data for the cell.
/// @param pointsCurvilinearCoordinates Evaluation points on the edge, in curvilinear coordinates.
/// @return A matrix of velocity edge basis-function values at the requested points.
Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                                   const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                   const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                                   const Eigen::MatrixXd &pointsCurvilinearCoordinates);

/// @brief Compute the (VEM) velocity stabilization matrix of the local space.
///
/// Returns the stabilization term complementing the consistency part of the VEM
/// velocity bilinear form; empty/zero for the Taylor–Hood FEM pair.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         Projection used to define the stabilization; defaults to the energy projection
/// \f$\Pi^\nabla\f$.
/// @return The local velocity stabilization matrix.
Eigen::MatrixXd VelocityStabilizationMatrix(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

/// @brief Coordinates of the velocity DOFs associated with a given edge.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_local_index       Local index of the edge within the cell.
/// @return The physical coordinates of the edge velocity DOFs (one per column).
Eigen::MatrixXd VelocityEdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                            const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                            const unsigned int edge_local_index);

/// @brief Number of local velocity degrees of freedom of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The velocity local-space dimension (number of local velocity DOFs).
unsigned int VelocitySize(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                          const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Compute the per-cell performance metrics of the local space.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The collected @ref Performance_Data (quadrature counts and VEM analysis).
Polydim::PDETools::LocalSpace_DF_PCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Export the velocity DOFs and solution fields to file for visualization.
///
/// Writes the right-hand side, the computed solution and the strongly-imposed
/// (Dirichlet) solution over the velocity DOFs to @p file_path, for post-processing.
/// The DOF layout and counts are provided per field.
///
/// @param geometry_utilities  GeDiM geometry helper.
/// @param mesh                Mesh data access object.
/// @param mesh_geometric_data Precomputed geometric data of the mesh.
/// @param mesh_dofs_info      Per-field mesh DOFs layout information.
/// @param dofs_data           Per-field DOF numbering/data.
/// @param count_dofs          Aggregated DOF counts across fields.
/// @param right_hand_side     Right-hand-side vector.
/// @param solution            Computed solution over the free DOFs.
/// @param solution_strongs    Values of the strongly-imposed (Dirichlet) DOFs.
/// @param file_path           Output file path.
void export_velocity_dofs(const Gedim::GeometryUtilities &geometry_utilities,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const std::vector<DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Gedim::IArray &right_hand_side,
                          const Gedim::IArray &solution,
                          const Gedim::IArray &solution_strongs,
                          const std::string &file_path);

} // namespace LocalSpace_DF_PCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
