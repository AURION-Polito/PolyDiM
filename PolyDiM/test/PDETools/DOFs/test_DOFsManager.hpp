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
//  struct Test_DOFsManager_mesh_connectivity_data final
//  {
//      inline unsigned int Cell0Ds_number() const
//      {
//          return mesh_data.Cell0DTotalNumber();
//      }
//      inline unsigned int Cell1Ds_number() const
//      {
//          return mesh_data.Cell1DTotalNumber();
//      }
//      inline unsigned int Cell2Ds_number() const
//      {
//          return mesh_data.Cell2DTotalNumber();
//      }
//      inline unsigned int Cell3Ds_number() const
//      {
//          return mesh_data.Cell3DTotalNumber();
//      }

//      inline unsigned int Cell0D_marker(const unsigned int cell0D_index) const
//      {
//          return mesh_data.Cell0DMarker(cell0D_index);
//      }
//      inline unsigned int Cell1D_marker(const unsigned int cell1D_index) const
//      {
//          return mesh_data.Cell1DMarker(cell1D_index);
//      }
//      inline unsigned int Cell2D_marker(const unsigned int cell2D_index) const
//      {
//          return mesh_data.Cell2DMarker(cell2D_index);
//      }
//      inline unsigned int Cell3D_marker(const unsigned int cell3D_index) const
//      {
//          return mesh_data.Cell3DMarker(cell3D_index);
//      }

//      inline std::array<unsigned int, 2> Cell1D_vertices(const unsigned int cell1D_index) const
//      {
//          return {mesh_data.Cell1DOrigin(cell1D_index), mesh_data.Cell1DEnd(cell1D_index)};
//      }

//      inline std::vector<unsigned int> Cell2D_vertices(const unsigned int cell2D_index) const
//      {
//          return mesh_data.Cell2DVertices(cell2D_index);
//      }

//      inline std::vector<unsigned int> Cell2D_edges(const unsigned int cell2D_index) const
//      {
//          return mesh_data.Cell2DEdges(cell2D_index);
//      }

//      inline std::vector<unsigned int> Cell3D_vertices(const unsigned int cell3D_index) const
//      {
//          return mesh_data.Cell3DVertices(cell3D_index);
//      }

//      inline std::vector<unsigned int> Cell3D_edges(const unsigned int cell3D_index) const
//      {
//          return mesh_data.Cell3DEdges(cell3D_index);
//      }

//      inline std::vector<unsigned int> Cell3D_faces(const unsigned int cell3D_index) const
//      {
//          return mesh_data.Cell3DFaces(cell3D_index);
//      }
//  };

TEST(TEST_DOFsManager, TEST_CreateDOFs_0D)
{
    Polydim::PDETools::DOFs::DOFsManager dofs_manager;

    PDETools::DOFs::DOFsManager::MeshDOFsInfo dofs_info;
    dofs_info.CellsBoundaryInfo[0].resize(1, {PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None, 0});
    dofs_info.CellsNumDOFs[0].resize(1, 2);

    const auto dofs_data = dofs_manager.CreateDOFs<0>(dofs_info);

    ASSERT_EQ(2, dofs_data.NumberDOFs);
    ASSERT_EQ(0, dofs_data.NumberStrongs);
    ASSERT_EQ(2, dofs_data.NumberInternalDOFs);
    ASSERT_EQ(0, dofs_data.NumberBoundaryDOFs);
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

    const auto dofs_data = dofs_manager.CreateDOFs<1>(dofs_info, Test_DOFsManager_mesh_connectivity_data_1D());

    ASSERT_EQ(3, dofs_data.NumberDOFs);
    ASSERT_EQ(1, dofs_data.NumberStrongs);
    ASSERT_EQ(2, dofs_data.NumberInternalDOFs);
    ASSERT_EQ(1, dofs_data.NumberBoundaryDOFs);
}

} // namespace UnitTesting
} // namespace Polydim

#endif
