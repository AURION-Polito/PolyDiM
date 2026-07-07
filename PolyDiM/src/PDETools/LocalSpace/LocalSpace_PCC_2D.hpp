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

#ifndef __LocalSpace_PCC_2D_H
#define __LocalSpace_PCC_2D_H

#include "DOFsManager.hpp"
#include "FEM_PCC_2D_LocalSpace.hpp"
#include "IArray.hpp"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_ZFEM_PCC_2D_LocalSpace.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_2D_Creator.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "ZFEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace PDETools
{

/// @brief Unified interface to the 2D primal conforming (PCC) local spaces.
///
/// This namespace provides a single, method-agnostic API for building and
/// evaluating the local (per-cell) discrete space of a scalar primal problem in 2D,
/// hiding the differences between the finite element (FEM), virtual element (VEM,
/// including its inertia-based and orthogonalized variants), and enriched FEM (ZFEM)
/// backends. The active backend is selected through MethodTypes; the free functions
/// dispatch to the corresponding implementation, so that assemblers can be written
/// once and reused across methods.
namespace LocalSpace_PCC_2D
{
/// @brief Discretization method for the 2D PCC local space.
enum struct MethodTypes
{
    FEM_PCC = 0,         ///< Finite Element Method.
    VEM_PCC = 1,         ///< Virtual Element Method (standard basis).
    VEM_PCC_Inertia = 2, ///< VEM with an inertia-based (principal-axes) monomial basis.
    VEM_PCC_Ortho = 3,   ///< VEM with an \f$L^2\f$-orthonormalized monomial basis.
    ZFEM_PCC = 4,        ///< Zipped Finite Element Method.
};

/// @brief Method-specific reference-element data for the 2D PCC local space.
///
/// Holds the selected method and polynomial order together with the reference-element
/// and local-space handles of the active backend. Only the members belonging to
/// @ref Method_Type are populated; the others remain unset.
class ReferenceElement_Data final
{
  public:
    Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes Method_Type; ///< Selected discretization method.
    unsigned int Order;                                            ///< Polynomial order of the space.

    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement> VEM_ReferenceElement;
    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement_Data;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_LocalSpace> VEM_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement> FEM_ReferenceElement;
    Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace> FEM_LocalSpace;

    std::unique_ptr<Polydim::ZFEM::PCC::I_ZFEM_PCC_2D_ReferenceElement> ZFEM_ReferenceElement;
    Polydim::ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data ZFEM_ReferenceElement_Data;
    std::unique_ptr<ZFEM::PCC::I_ZFEM_PCC_2D_LocalSpace> ZFEM_LocalSpace;
};

/// @brief Per-cell local-space data for the 2D PCC local space.
///
/// Stores the polygon geometry and the computed local-space data for the active
/// backend; only the members of the selected method are populated.
class LocalSpace_Data final
{
  public:
    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data VEM_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data FEM_LocalSpace_Data;

    Polydim::ZFEM::PCC::ZFEM_PCC_2D_Polygon_Geometry ZFEM_Geometry;
    Polydim::ZFEM::PCC::ZFEM_PCC_2D_LocalSpace_Data ZFEM_LocalSpace_Data;
};

/// @brief Performance metrics collected on a single 2D cell.
class Performance_Data final
{
  public:
    class Cell2D_Performance final
    {
      public:
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data vem_analysis_data;
        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis_Data zfem_analysis_data;
    };

    Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data::Cell2D_Performance performance_data;
};

/// @brief Create the reference element for the chosen method and order.
///
/// Factory that instantiates the backend-specific reference element and populates the
/// corresponding members of the returned @ref ReferenceElement_Data.
///
/// @param method_type  Discretization method to use.
/// @param method_order Polynomial order of the space.
/// @return The initialized reference-element data.
Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes &method_type,
                                                                                   const unsigned int method_order);

/// @brief Geometric-data configuration required by the selected method.
///
/// Returns the GeDiM mesh geometric-data configuration (which geometric quantities must
/// be precomputed on each cell) needed by the active backend.
///
/// @param reference_element_data Reference-element data describing the method.
/// @return The mesh geometric-data configuration.
/// @note The spelling of this function name (@c MeshGeometricDataConfigiguration) is preserved from the source.
Gedim::MeshUtilities::MeshGeometricData2DConfig MeshGeometricDataConfigiguration(const ReferenceElement_Data &reference_element_data);

/// @brief Build the degrees-of-freedom layout over the mesh.
///
/// Determines how the DOFs of the chosen space are distributed over mesh entities
/// (vertices, edges, cells) and marks boundary DOFs according to @p boundary_info.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param mesh                   Mesh data access object.
/// @param boundary_info          Per-marker boundary information (Dirichlet/Neumann, etc.).
/// @return The mesh DOFs layout information.
Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo SetMeshDOFsInfo(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> &boundary_info);

/// @brief Build the local space on a given cell.
///
/// Computes the backend-specific local-space data on cell @p cell2D_index from the
/// precomputed mesh geometric data, returning the @ref LocalSpace_Data used by all
/// subsequent evaluation routines.
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
                                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data);

/// @brief Evaluate the basis functions at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM/ZFEM); defaults to \f$\Pi^0_{k-1}\f$.
/// @return A (quadrature points \f$\times\f$ local DOFs) matrix of basis-function values.
Eigen::MatrixXd BasisFunctionsValues(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                     const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

/// @brief Evaluate the basis functions at arbitrary points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param points                 Evaluation points (one per column).
/// @param projectionType         VEM projection operator to apply (ignored by FEM/ZFEM); defaults to \f$\Pi^0_{k-1}\f$.
/// @return A (points \f$\times\f$ local DOFs) matrix of basis-function values.
Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                     const LocalSpace_Data &local_space_data,
                                     const Eigen::MatrixXd &points,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

/// @brief Evaluate the basis functions on a given edge.
///
/// @param edge_local_index           Local index of the edge within the cell.
/// @param reference_element_data     Reference-element data describing the method.
/// @param local_space_data           Local-space data for the cell.
/// @param pointsCurvilinearCoordinates Evaluation points on the edge, in curvilinear coordinates.
/// @return A matrix of edge basis-function values at the requested points.
Eigen::MatrixXd BasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &pointsCurvilinearCoordinates);

/// @brief Evaluate the basis-function derivatives at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM/ZFEM); defaults to the derivative
/// projection \f$\Pi^0_{k-1}\f$.
/// @return One (quadrature points \f$\times\f$ local DOFs) matrix per spatial direction.
std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

/// @brief Evaluate the basis-function derivatives at arbitrary points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param points                 Evaluation points (one per column).
/// @param projectionType         VEM projection operator to apply (ignored by FEM/ZFEM); defaults to the derivative
/// projection \f$\Pi^0_{k-1}\f$.
/// @return One (points \f$\times\f$ local DOFs) matrix per spatial direction.
std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
    const Eigen::MatrixXd &points,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

/// @brief Evaluate the basis-function Laplacians at the internal quadrature points.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         VEM projection operator to apply (ignored by FEM/ZFEM); defaults to the derivative
/// projection \f$\Pi^0_{k-1}\f$.
/// @return A (quadrature points \f$\times\f$ local DOFs) matrix of Laplacian values.
Eigen::MatrixXd BasisFunctionsLaplacianValues(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

/// @brief Compute the (VEM) stabilization matrix of the local space.
///
/// Returns the stabilization term that complements the consistency part of the VEM
/// bilinear form; for FEM/ZFEM, where no stabilization is required, the term is empty
/// or zero as appropriate.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param projectionType         Projection used to define the stabilization; defaults to the energy projection
/// \f$\Pi^\nabla\f$.
/// @return The local stabilization matrix.
Eigen::MatrixXd StabilizationMatrix(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::PiNabla);

/// @brief Coordinates of the DOFs associated with a given edge.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @param edge_local_index       Local index of the edge within the cell.
/// @return The physical coordinates of the edge DOFs (one per column).
Eigen::MatrixXd EdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index);

/// @brief Coordinates (and weights) of the internal DOFs of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return Quadrature-like data holding the internal-DOF coordinates and weights.
Gedim::Quadrature::QuadratureData InternalDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                          const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Compute the internal DOF values of a function.
///
/// Given the function values sampled at the internal-DOF coordinates, returns the
/// corresponding internal degrees of freedom (e.g. the interior moments for VEM).
///
/// @param reference_element_data    Reference-element data describing the method.
/// @param local_space_data          Local-space data for the cell.
/// @param values_at_dofs            Function values at the internal-DOF coordinates.
/// @param internal_dofs_coordinates Internal-DOF coordinates/weights from InternalDofsCoordinates().
/// @return The vector of internal degrees of freedom.
Eigen::VectorXd InternalDofs(const ReferenceElement_Data &reference_element_data,
                             const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                             const Eigen::VectorXd &values_at_dofs,
                             const Gedim::Quadrature::QuadratureData &internal_dofs_coordinates);

/// @brief Internal quadrature rule of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The internal quadrature points and weights used for volume integration.
Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Number of local degrees of freedom of the cell.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The local-space dimension (number of local DOFs).
unsigned int Size(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                  const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Compute the per-cell performance metrics of the local space.
///
/// @param reference_element_data Reference-element data describing the method.
/// @param local_space_data       Local-space data for the cell.
/// @return The collected @ref Performance_Data (quadrature counts and method-specific analysis).
Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

/// @brief Export the DOFs and solution fields to file for visualization.
///
/// Writes the right-hand side, the computed solution and the strongly-imposed
/// (Dirichlet) solution over the mesh DOFs to @p file_path, for post-processing.
///
/// @param geometry_utilities  GeDiM geometry helper.
/// @param mesh                Mesh data access object.
/// @param mesh_geometric_data Precomputed geometric data of the mesh.
/// @param mesh_dofs_info      Mesh DOFs layout information.
/// @param dofs_data           DOF numbering/data.
/// @param right_hand_side     Right-hand-side vector.
/// @param solution            Computed solution over the free DOFs.
/// @param solution_strongs    Values of the strongly-imposed (Dirichlet) DOFs.
/// @param file_path           Output file path.
void export_dofs(const Gedim::GeometryUtilities &geometry_utilities,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Gedim::IArray &right_hand_side,
                 const Gedim::IArray &solution,
                 const Gedim::IArray &solution_strongs,
                 const std::string &file_path);

} // namespace LocalSpace_PCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
