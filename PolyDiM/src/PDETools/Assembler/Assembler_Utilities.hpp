#ifndef __PDETOOLS_ASSEMBLER_AssemblerUtilities_HPP
#define __PDETOOLS_ASSEMBLER_AssemblerUtilities_HPP

#include "DOFsManager.hpp"
#include <vector>

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
struct local_matrix_to_global_matrix_dofs_data final
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data;
    std::vector<unsigned int> local_offsets;
    std::vector<unsigned int> global_offsets_DOFs;
    std::vector<unsigned int> global_offsets_Strongs;
};

template <unsigned int dimension,
          typename local_lhs_type,
          typename local_rhs_type,
          typename global_lhs_type,
          typename global_rhs_type>
void assemble_local_matrix_to_global_matrix(const unsigned int cell_index,
                                            const local_matrix_to_global_matrix_dofs_data &test_functions_dofs_data,
                                            const local_matrix_to_global_matrix_dofs_data &trial_functions_dofs_data,
                                            const local_lhs_type &local_lhs,
                                            const local_rhs_type &local_rhs,
                                            global_lhs_type &global_lhs_DOFs,
                                            global_lhs_type &global_lhs_Strongs,
                                            global_rhs_type &global_rhs)
{
    for (unsigned test_f = 0; test_f < test_functions_dofs_data.dofs_data.size(); ++test_f)
    {
        const auto &test_dofs_data = test_functions_dofs_data.dofs_data.at(test_f).get();
        const auto &test_global_dofs = test_dofs_data.CellsGlobalDOFs.at(dimension).at(cell_index);
        const auto test_global_offset_DOFs = test_functions_dofs_data.global_offsets_DOFs[test_f];
        const auto test_local_offset = test_functions_dofs_data.local_offsets[test_f];

        for (unsigned int test_loc_i = 0; test_loc_i < test_global_dofs.size(); ++test_loc_i)
        {
            const auto global_dof_i = test_global_dofs.at(test_loc_i);
            const auto local_dof_i = test_dofs_data.CellsDOFs.at(global_dof_i.Dimension)
                                         .at(global_dof_i.CellIndex)
                                         .at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                continue;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }

            const unsigned int global_index_i = local_dof_i.Global_Index + test_global_offset_DOFs;

            global_rhs.AddValue(global_index_i, local_rhs[test_loc_i + test_local_offset]);

            for (unsigned int trial_f = 0; trial_f < trial_functions_dofs_data.dofs_data.size(); ++trial_f)
            {
                const auto &trial_dofs_data = trial_functions_dofs_data.dofs_data.at(trial_f).get();
                const auto &trial_global_dofs = trial_dofs_data.CellsGlobalDOFs.at(dimension).at(cell_index);
                const auto trial_global_offset_DOFs = trial_functions_dofs_data.global_offsets_DOFs[trial_f];
                const auto trial_local_offset = trial_functions_dofs_data.local_offsets[trial_f];
                const auto trial_global_offset_Strongs = trial_functions_dofs_data.global_offsets_Strongs[trial_f];

                for (unsigned int trial_loc_j = 0; trial_loc_j < trial_global_dofs.size(); ++trial_loc_j)
                {
                    const auto &global_dof_j = trial_global_dofs.at(trial_loc_j);
                    const auto &local_dof_j = trial_dofs_data.CellsDOFs.at(global_dof_j.Dimension)
                                                  .at(global_dof_j.CellIndex)
                                                  .at(global_dof_j.DOFIndex);

                    const unsigned int global_index_j = local_dof_j.Global_Index;
                    const double loc_A_element =
                        local_lhs(test_loc_i + test_local_offset, trial_loc_j + trial_local_offset);

                    switch (local_dof_j.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        global_lhs_Strongs.Triplet(
                            global_index_i, global_index_j + trial_global_offset_Strongs, loc_A_element);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        global_lhs_DOFs.Triplet(
                            global_index_i, global_index_j + trial_global_offset_DOFs, loc_A_element);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }
        }
    }
}
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
