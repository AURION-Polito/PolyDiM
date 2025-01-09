#ifndef __program_utilities_H
#define __program_utilities_H

#include "DOFsManager.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"
#include <typeindex>
#include <unordered_map>

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
namespace program_utilities
{
unique_ptr<Polydim::examples::Stokes_DF_PCC_3D::test::I_Test> create_test(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config);

void create_domain_mesh(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &domain,
                        Gedim::MeshMatricesDAO &mesh);

Gedim::MeshUtilities::MeshGeometricData3D create_domain_mesh_geometric_properties(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh);

void export_solution(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder);

void export_velocity_dofs(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &vem_velocity_reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder);

void export_performance(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                        const Assembler::VEM_Performance_Result &performance_data,
                        const std::string &exportFolder);

void export_discrepancy_errors(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Stokes_DF_PCC_3D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
                               const std::string &exportSolutionFolder,
                               const std::string &exportVtuFolder);

} // namespace program_utilities
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
