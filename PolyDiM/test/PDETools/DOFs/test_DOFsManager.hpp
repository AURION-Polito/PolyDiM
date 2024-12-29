#ifndef __TEST_CreateDOFs_H
#define __TEST_CreateDOFs_H

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

    const auto dofs_data = dofs_manager.CreateDOFs<0>(dofs_info);

    ASSERT_EQ(1, dofs_data.NumberDOFs);
    ASSERT_EQ(0, dofs_data.NumberStrongs);
    ASSERT_EQ(1, dofs_data.NumberInternalDOFs);
    ASSERT_EQ(0, dofs_data.NumberBoundaryDOFs);
}

} // namespace UnitTesting
} // namespace Polydim

#endif
