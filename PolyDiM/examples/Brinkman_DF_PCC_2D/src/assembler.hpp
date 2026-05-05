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

    struct Performance_Data final
    {
        std::vector<Polydim::PDETools::LocalSpace_DF_PCC_2D::Performance_Data> Cell2DsPerformance;
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
    void ComputeStrongTerm(const unsigned int cell2D_index,
                           const Gedim::MeshMatricesDAO &mesh,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                           const std::vector<size_t> &offsetStrongs,
                           const PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                           const PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                           const test::I_Test &test,
                           Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                         const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                         const PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                         const PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                         const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                         Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

    std::map<unsigned int, double> ComputeFlux(const Program_configuration &config,
                                               const Gedim::MeshMatricesDAO &mesh,
                                               const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                               const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                               const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                               const PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                               const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                                               const Stokes_DF_PCC_2D_Problem_Data &assembler_data) const;

  public:
    Stokes_DF_PCC_2D_Problem_Data Assemble(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                           const Gedim::MeshMatricesDAO &mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                           const Polydim::PDETools::Assembler_Utilities::count_dofs_data &counts_dofs,
                                           const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                           const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const;

    Performance_Data ComputeMethodPerformance(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                              const Gedim::MeshMatricesDAO &mesh,
                                              const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                              const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                         const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                         const Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const;

    Assembler::DiscrepancyErrors_Data ComputeDiscrepancyErrors(
        const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const std::vector<PDETools::DOFs::DOFsManager::DOFsData> &reduced_dofs_data,
        const Polydim::PDETools::Assembler_Utilities::count_dofs_data &reduced_count_dofs,
        const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reduced_reference_element_data,
        const Stokes_DF_PCC_2D_Problem_Data &reduced_assembler_data,
        const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const;
};
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
