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
namespace Parabolic_PCC_2D
{
namespace program_utilities
{

std::unique_ptr<Polydim::examples::Parabolic_PCC_2D::test::I_Test> create_test(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config);

void create_domain_mesh(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                        Gedim::MeshMatricesDAO &mesh);

Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh);

void export_solution(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Gedim::Eigen_SparseArray<> &A,
                     const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const unsigned int time_index,
                     const double &time_value,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder);

void export_dofs(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Problem_Data &assembler_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                 const unsigned int time_index,
                 const double &time_value,
                 const std::string &exportVtuFolder);

void export_performance(const Program_configuration &config,
                        const Assembler::Performance_Data &performance_data,
                        const std::string &exportFolder);

std::vector<double> create_time_steps(const Program_configuration &config, const std::array<double, 2> &time_domain);

} // namespace program_utilities
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
