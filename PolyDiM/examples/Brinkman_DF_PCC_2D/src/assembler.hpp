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

#ifndef __assembler_H
#define __assembler_H

#include "Assembler_Utilities.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

#include "DOFsManager.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"

#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
class Assembler final
{
  public:
    struct Stokes_DF_PCC_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionDirichlet;
    };

    struct VEM_Performance_Result final
    {
        struct Cell2D_Performance final
        {
            unsigned int NumBoundaryQuadraturePoints = 0;
            unsigned int NumInternalQuadraturePoints = 0;
            double maxPiNablaConditioning;    ///< conditioning of piNabla
            double maxPi0kConditioning;       ///< conditioning of piNabla
            double maxErrorPiNabla;           ///< |piNabla * Dofs - I|
            double maxErrorPi0k;              ///< |pi0k * Dofs - I|
            double ErrorStabilization = -1.0; ///< |S * Dofs|
            double maxErrorHCD;               ///< |H - CD|
            double maxErrorGBD;               ///< |G - BD|
        };

        std::vector<Cell2D_Performance> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        std::array<Eigen::VectorXd, 3> cell0Ds_numeric_velocity;
        std::array<Eigen::VectorXd, 3> cell0Ds_exact_velocity;

        Eigen::VectorXd cell2Ds_discrepancy_error_H1_velocity;
        Eigen::VectorXd cell2Ds_error_H1_velocity;
        Eigen::VectorXd cell2Ds_norm_H1_velocity;
        double discrepancy_error_H1_velocity;
        double error_H1_velocity;
        double norm_H1_velocity;

        Eigen::MatrixXd repeated_vertices_coordinates;
        std::vector<std::vector<unsigned int>> repeated_connectivity;
        Eigen::VectorXd cell0Ds_numeric_pressure;
        Eigen::VectorXd cell0Ds_exact_pressure;

        Eigen::VectorXd cell2Ds_discrepancy_error_L2_pressure;
        Eigen::VectorXd cell2Ds_error_L2_pressure;
        Eigen::VectorXd cell2Ds_norm_L2_pressure;
        double discrepancy_error_L2_pressure;
        double error_L2_pressure;
        double norm_L2_pressure;

        double mesh_size;

        double residual_norm;
        std::map<unsigned int, double> flux;

        Eigen::VectorXd inverse_diffusion_coeff_values;
        Eigen::VectorXd viscosity_values;
    };

    struct DiscrepancyErrors_Data final
    {
        Eigen::VectorXd cell2Ds_discrepancy_error_L2_pressure;
        Eigen::VectorXd cell2Ds_discrepancy_error_H1_velocity;
        Eigen::VectorXd cell2Ds_full_norm_L2_pressure;
        Eigen::VectorXd cell2Ds_full_norm_H1_velocity;
        double discrepancy_error_L2_pressure;
        double discrepancy_error_H1_velocity;
        double full_norm_L2_pressure;
        double full_norm_H1_velocity;

        double residual_norm;
        double reduced_residual_norm;
        double pressure_dofs_ratio;
        double velocity_dofs_ratio;
    };

  private:
    void ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                           const std::vector<size_t> &offsetStrongs,
                           const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                           const test::I_Test &test,
                           Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                         const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                         const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_local_space,
                         const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                         Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

    std::map<unsigned int, double> ComputeFlux(const Program_configuration &config,
                                               const Gedim::MeshMatricesDAO &mesh,
                                               const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                               const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                               const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                               const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                               const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_local_space,
                                               const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                                               const Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

  public:
    Stokes_DF_PCC_2D_Problem_Data Assemble(
        const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
        const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
        const Polydim::PDETools::Assembler_Utilities::count_dofs_data &counts_dofs,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
        const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const;

    VEM_Performance_Result ComputeVemPerformance(
        const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                         const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                         const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
                                         const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
                                         const Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const;

    Assembler::DiscrepancyErrors_Data ComputeDiscrepancyErrors(
        const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &full_dofs_data,
        const Polydim::PDETools::Assembler_Utilities::count_dofs_data &full_count_dofs,
        const std::vector<PDETools::DOFs::DOFsManager::DOFsData> &reduced_dofs_data,
        const Polydim::PDETools::Assembler_Utilities::count_dofs_data &reduced_count_dofs,
        const VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &full_velocity_reference_element_data,
        const VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &full_pressure_reference_element_data,
        const VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reduced_velocity_reference_element_data,
        const VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reduced_pressure_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_full_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_full_pressure_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_reduced_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_reduced_pressure_local_space,
        const Stokes_DF_PCC_2D_Problem_Data &full_assembler_data,
        const Stokes_DF_PCC_2D_Problem_Data &reduced_assembler_data) const;
};
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
