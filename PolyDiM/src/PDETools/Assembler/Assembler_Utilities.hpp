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
    std::vector<size_t> local_offsets;
    std::vector<size_t> global_offsets_DOFs;
    std::vector<size_t> global_offsets_Strongs;
};

struct count_dofs_data final
{
    unsigned int num_total_boundary_dofs;
    unsigned int num_total_dofs;
    unsigned int num_total_strong;
    std::vector<size_t> offsets_DOFs;
    std::vector<size_t> offsets_Strongs;
};

struct local_count_dofs_data final
{
    unsigned int num_total_dofs;
    std::vector<size_t> offsets_DOFs;
};

inline Polydim::PDETools::Assembler_Utilities::count_dofs_data count_dofs(const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data)
{
    Polydim::PDETools::Assembler_Utilities::count_dofs_data data;
    data.num_total_dofs = dofs_data[0].NumberDOFs;
    data.num_total_strong = dofs_data[0].NumberStrongs;
    data.num_total_boundary_dofs = dofs_data[0].NumberBoundaryDOFs;

    const unsigned int numDOFHandler = dofs_data.size();
    data.offsets_DOFs.resize(numDOFHandler);
    data.offsets_Strongs.resize(numDOFHandler);
    data.offsets_DOFs[0] = 0;
    data.offsets_Strongs[0] = 0;

    for (unsigned int i = 1; i < numDOFHandler; i++)
    {
        data.num_total_dofs += dofs_data[i].NumberDOFs;
        data.num_total_strong += dofs_data[i].NumberStrongs;
        data.num_total_boundary_dofs += dofs_data[i].NumberBoundaryDOFs;

        data.offsets_DOFs[i] = data.offsets_DOFs[i - 1] + dofs_data[i - 1].NumberDOFs;
        data.offsets_Strongs[i] = data.offsets_Strongs[i - 1] + dofs_data[i - 1].NumberStrongs;
    }

    return data;
}
// ***************************************************************************
template <unsigned int dimension>
inline Polydim::PDETools::Assembler_Utilities::local_count_dofs_data local_count_dofs(
    const unsigned int cell_index,
    const std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> &dofs_data)
{
    Polydim::PDETools::Assembler_Utilities::local_count_dofs_data data;
    data.num_total_dofs = dofs_data[0].get().CellsGlobalDOFs[dimension].at(cell_index).size();

    const unsigned int numDOFHandler = dofs_data.size();
    data.offsets_DOFs.resize(numDOFHandler);
    data.offsets_DOFs[0] = 0;

    for (unsigned int i = 1; i < numDOFHandler; i++)
    {
        data.num_total_dofs += dofs_data[i].get().CellsGlobalDOFs[dimension].at(cell_index).size();
        data.offsets_DOFs[i] =
            data.offsets_DOFs[i - 1] + dofs_data[i - 1].get().CellsGlobalDOFs[dimension].at(cell_index).size();
    }

    return data;
}
// ***************************************************************************
template <unsigned int dimension>
inline Polydim::PDETools::Assembler_Utilities::local_count_dofs_data
local_count_dofs(const unsigned int cell_index, const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data)
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(1);
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));

    return local_count_dofs<dimension>(cell_index, dofs_data_ref);
}
// ***************************************************************************
template <unsigned int dimension>
inline Polydim::PDETools::Assembler_Utilities::local_count_dofs_data
local_count_dofs(const unsigned int cell_index, const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data)
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(dofs_data.size());
    for (unsigned int c = 0; c < dofs_data.size(); c++)
        dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data.at(c))));

    return local_count_dofs<dimension>(cell_index, dofs_data_ref);
}
// ***************************************************************************
template <unsigned int dimension, typename global_solution_type>
inline Eigen::VectorXd global_solution_to_local_solution(
    const unsigned int cell_index,
    const std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> &dofs_data,
    const unsigned int &num_local_DOFs,
    const std::vector<size_t> &local_offsets_DOFs,
    const std::vector<size_t> &global_offsets_DOFs,
    const std::vector<size_t> &global_offsets_Strongs,
    const global_solution_type &global_solution_DOFs,
    const global_solution_type &global_solution_Strongs)
{
    Eigen::VectorXd local_solution_dofs = Eigen::VectorXd::Zero(num_local_DOFs);

    const unsigned int numDOFHandler = dofs_data.size();
    for (unsigned int h = 0; h < numDOFHandler; h++)
    {
        const auto &dofs = dofs_data[h].get();
        const auto &global_dof = dofs.CellsGlobalDOFs[dimension].at(cell_index);
        for (unsigned int loc_i = 0; loc_i < global_dof.size(); ++loc_i)
        {
            const auto &global_dof_i = global_dof.at(loc_i);
            const auto &local_dof_i =
                dofs.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                local_solution_dofs[loc_i + local_offsets_DOFs[h]] =
                    global_solution_Strongs.GetValue(local_dof_i.Global_Index + global_offsets_Strongs[h]);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                local_solution_dofs[loc_i + local_offsets_DOFs[h]] =
                    global_solution_DOFs.GetValue(local_dof_i.Global_Index + global_offsets_DOFs[h]);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    return local_solution_dofs;
}
// ***************************************************************************
template <unsigned int dimension, typename global_solution_type>
inline Eigen::VectorXd global_solution_to_local_solution(const unsigned int cell_index,
                                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                         const unsigned int &num_local_DOFs,
                                                         const std::vector<size_t> &local_offsets_DOFs,
                                                         const std::vector<size_t> &global_offsets_DOFs,
                                                         const std::vector<size_t> &global_offsets_Strongs,
                                                         const global_solution_type &global_solution_DOFs,
                                                         const global_solution_type &global_solution_Strongs)
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(1);
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));

    return global_solution_to_local_solution<dimension, global_solution_type>(cell_index,
                                                                              dofs_data_ref,
                                                                              num_local_DOFs,
                                                                              local_offsets_DOFs,
                                                                              global_offsets_DOFs,
                                                                              global_offsets_Strongs,
                                                                              global_solution_DOFs,
                                                                              global_solution_Strongs);
}
// ***************************************************************************
template <unsigned int dimension, typename global_solution_type>
inline Eigen::VectorXd global_solution_to_local_solution(const unsigned int cell_index,
                                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                                         const unsigned int &num_local_DOFs,
                                                         const std::vector<size_t> &local_offsets_DOFs,
                                                         const std::vector<size_t> &global_offsets_DOFs,
                                                         const std::vector<size_t> &global_offsets_Strongs,
                                                         const global_solution_type &global_solution_DOFs,
                                                         const global_solution_type &global_solution_Strongs)
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(dofs_data.size());
    for (unsigned int c = 0; c < dofs_data.size(); c++)
        dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data.at(c))));

    return global_solution_to_local_solution<dimension, global_solution_type>(cell_index,
                                                                              dofs_data_ref,
                                                                              num_local_DOFs,
                                                                              local_offsets_DOFs,
                                                                              global_offsets_DOFs,
                                                                              global_offsets_Strongs,
                                                                              global_solution_DOFs,
                                                                              global_solution_Strongs);
}
// ***************************************************************************
template <unsigned int dimension, typename local_lhs_type, typename local_rhs_type, typename global_lhs_type, typename global_rhs_type>
void assemble_local_matrix_to_global_matrix(const unsigned int cell_index,
                                            const local_matrix_to_global_matrix_dofs_data &test_functions_dofs_data,
                                            const local_matrix_to_global_matrix_dofs_data &trial_functions_dofs_data,
                                            const local_lhs_type &local_lhs,
                                            const local_rhs_type &local_rhs,
                                            global_lhs_type &global_lhs_DOFs,
                                            global_lhs_type &global_lhs_Strongs,
                                            global_rhs_type &global_rhs)
{
    for (size_t test_f = 0; test_f < test_functions_dofs_data.dofs_data.size(); ++test_f)
    {
        const auto &test_dofs_data = test_functions_dofs_data.dofs_data.at(test_f).get();
        const auto &test_global_dofs = test_dofs_data.CellsGlobalDOFs.at(dimension).at(cell_index);
        const auto test_global_offset_DOFs = test_functions_dofs_data.global_offsets_DOFs[test_f];
        const auto test_local_offset = test_functions_dofs_data.local_offsets[test_f];

        for (size_t test_loc_i = 0; test_loc_i < test_global_dofs.size(); ++test_loc_i)
        {
            const auto global_dof_i = test_global_dofs.at(test_loc_i);
            const auto local_dof_i =
                test_dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

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
                    const auto &local_dof_j =
                        trial_dofs_data.CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

                    const unsigned int global_index_j = local_dof_j.Global_Index;
                    const double loc_A_element = local_lhs(test_loc_i + test_local_offset, trial_loc_j + trial_local_offset);

                    switch (local_dof_j.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        global_lhs_Strongs.Triplet(global_index_i, global_index_j + trial_global_offset_Strongs, loc_A_element);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        global_lhs_DOFs.Triplet(global_index_i, global_index_j + trial_global_offset_DOFs, loc_A_element);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }
        }
    }
}
// ***************************************************************************
template <unsigned int dimension, typename local_lhs_type, typename global_lhs_type>
void assemble_local_matrix_to_global_matrix(const unsigned int cell_index,
                                            const local_matrix_to_global_matrix_dofs_data &test_functions_dofs_data,
                                            const local_matrix_to_global_matrix_dofs_data &trial_functions_dofs_data,
                                            const local_lhs_type &local_lhs,
                                            global_lhs_type &global_lhs_DOFs,
                                            global_lhs_type &global_lhs_Strongs)
{
    for (size_t test_f = 0; test_f < test_functions_dofs_data.dofs_data.size(); ++test_f)
    {
        const auto &test_dofs_data = test_functions_dofs_data.dofs_data.at(test_f).get();
        const auto &test_global_dofs = test_dofs_data.CellsGlobalDOFs.at(dimension).at(cell_index);
        const auto test_global_offset_DOFs = test_functions_dofs_data.global_offsets_DOFs[test_f];
        const auto test_local_offset = test_functions_dofs_data.local_offsets[test_f];

        for (size_t test_loc_i = 0; test_loc_i < test_global_dofs.size(); ++test_loc_i)
        {
            const auto global_dof_i = test_global_dofs.at(test_loc_i);
            const auto local_dof_i =
                test_dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

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
                    const auto &local_dof_j =
                        trial_dofs_data.CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

                    const unsigned int global_index_j = local_dof_j.Global_Index;
                    const double loc_A_element = local_lhs(test_loc_i + test_local_offset, trial_loc_j + trial_local_offset);

                    switch (local_dof_j.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        global_lhs_Strongs.Triplet(global_index_i, global_index_j + trial_global_offset_Strongs, loc_A_element);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        global_lhs_DOFs.Triplet(global_index_i, global_index_j + trial_global_offset_DOFs, loc_A_element);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }
        }
    }
}
// ***************************************************************************
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
