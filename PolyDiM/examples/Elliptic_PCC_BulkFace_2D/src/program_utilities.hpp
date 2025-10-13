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
namespace Elliptic_PCC_BulkFace_2D
{
namespace program_utilities
{

std::unique_ptr<Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test> create_test(
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config);

void create_domain_mesh_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                           const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                           Gedim::MeshMatricesDAO &mesh);

std::vector<double> create_time_steps(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                      const std::array<double, 2> &time_domain);

void export_solution(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                     const double &value_time,
                     const Gedim::MeshMatricesDAO &mesh_2D,
                     const Gedim::MeshMatricesDAO &mesh_1D,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder);

void export_performance_2D(const Program_configuration &config,
                           const Assembler::Performance_Data_2D &performance_data,
                           const std::string &exportFolder);

void export_solution_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                        const double &value_time,
                        const Gedim::MeshMatricesDAO &mesh,
                        const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                        const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data_1D &post_process_data,
                        const std::string &exportSolutionFolder,
                        const std::string &exportVtuFolder);

void export_solution_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                        const double &value_time,
                        const Gedim::MeshMatricesDAO &mesh,
                        const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                        const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data_2D &post_process_data,
                        const std::string &exportSolutionFolder,
                        const std::string &exportVtuFolder);

} // namespace program_utilities
} // namespace Elliptic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim

#endif
