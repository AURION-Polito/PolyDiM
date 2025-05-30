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

#ifndef __program_utilities_H
#define __program_utilities_H

#include "DOFsManager.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
namespace program_utilities
{
std::unique_ptr<Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test> create_test(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config);

void create_domain_mesh(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                        Gedim::MeshMatricesDAO &mesh);

Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh);

void export_solution(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder);

void export_velocity_dofs(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &vem_velocity_reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                          const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder);

void export_discrepancy_errors(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
                               const std::string &exportSolutionFolder,
                               const std::string &exportVtuFolder);

} // namespace program_utilities
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
