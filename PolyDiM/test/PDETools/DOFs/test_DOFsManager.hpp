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

#ifndef __TEST_DOFsManager_H
#define __TEST_DOFsManager_H

#include "DOFsManager.hpp"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace Polydim
{
namespace UnitTesting
{

TEST(TEST_DOFsManager, TEST_CreateDOFs_0D)
{
    Polydim::PDETools::DOFs::DOFsManager dofs_manager;

    PDETools::DOFs::DOFsManager::MeshDOFsInfo dofs_info;
    dofs_info.CellsBoundaryInfo[0].resize(1, {PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None, 0});
    dofs_info.CellsNumDOFs[0].resize(1, 2);

    const auto dofs_data = dofs_manager.CreateDOFs_0D(dofs_info);

    ASSERT_EQ(2, dofs_data.NumberDOFs);
    ASSERT_EQ(0, dofs_data.NumberStrongs);
    ASSERT_EQ(2, dofs_data.NumberInternalDOFs);
    ASSERT_EQ(0, dofs_data.NumberBoundaryDOFs);

    const auto cells_dofs_indices = dofs_manager.ComputeCellsDOFsIndices(dofs_data,
                                                                         0);

    ASSERT_EQ(std::vector<std::vector<unsigned int>>({ { 0u, 1u } }), cells_dofs_indices.Cells_DOFs_LocalIndex);
    ASSERT_EQ(std::vector<std::vector<unsigned int>>({ { 0u, 1u } }), cells_dofs_indices.Cells_DOFs_GlobalIndex);
    ASSERT_EQ(std::vector<std::vector<unsigned int>>({ { } }), cells_dofs_indices.Cells_Strongs_LocalIndex);
    ASSERT_EQ(std::vector<std::vector<unsigned int>>({ { } }), cells_dofs_indices.Cells_Strongs_GlobalIndex);
}

TEST(TEST_DOFsManager, TEST_CreateDOFs_1D)
{
    struct Test_DOFsManager_mesh_connectivity_data_1D final
    {
        inline unsigned int Cell0Ds_number() const
        {
            return 2;
        }
        inline unsigned int Cell1Ds_number() const
        {
            return 1;
        }

        inline unsigned int Cell0D_marker(const unsigned int cell0D_index) const
        {
            throw std::runtime_error("Not implemented");
        }
        inline unsigned int Cell1D_marker(const unsigned int cell1D_index) const
        {
            throw std::runtime_error("Not implemented");
        }

        inline std::array<unsigned int, 2> Cell1D_vertices(const unsigned int cell1D_index) const
        {
            return {0, 1};
        }
    };

    Polydim::PDETools::DOFs::DOFsManager dofs_manager;

    PDETools::DOFs::DOFsManager::MeshDOFsInfo dofs_info;
    dofs_info.CellsBoundaryInfo[0].resize(2);
    dofs_info.CellsBoundaryInfo[0][0] = {PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong, 2};
    dofs_info.CellsBoundaryInfo[0][1] = {PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak, 1};
    dofs_info.CellsNumDOFs[0].resize(2, 1);
    dofs_info.CellsBoundaryInfo[1].resize(1, {PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None, 0});
    dofs_info.CellsNumDOFs[1].resize(1, 2);

    const auto dofs_data = dofs_manager.CreateDOFs_1D(dofs_info, Test_DOFsManager_mesh_connectivity_data_1D());

    ASSERT_EQ(3, dofs_data.NumberDOFs);
    ASSERT_EQ(1, dofs_data.NumberStrongs);
    ASSERT_EQ(2, dofs_data.NumberInternalDOFs);
    ASSERT_EQ(1, dofs_data.NumberBoundaryDOFs);

    const auto cells_dofs_indices = dofs_manager.ComputeCellsDOFsIndices(dofs_data,
                                                                         1);
}

} // namespace UnitTesting
} // namespace Polydim

#endif
