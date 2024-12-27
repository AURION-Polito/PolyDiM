#ifndef __program_utilities_H
#define __program_utilities_H

#include "DOFsManager.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
namespace program_utilities
{

std::unique_ptr<Polydim::examples::Elliptic_PCC_3D::test::I_Test> create_test(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config);

void create_domain_mesh(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &domain,
                        Gedim::MeshMatricesDAO &mesh);

Gedim::MeshUtilities::MeshGeometricData3D create_domain_mesh_geometric_properties(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh);

void export_solution(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Polydim::examples::Elliptic_PCC_3D::Assembler::Elliptic_PCC_3D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder);

} // namespace program_utilities
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
