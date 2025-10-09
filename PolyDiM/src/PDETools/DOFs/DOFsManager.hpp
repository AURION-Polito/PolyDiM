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

#ifndef __PDETOOLS_DOFS_DOFsManager_HPP
#define __PDETOOLS_DOFS_DOFsManager_HPP

#include "Eigen/Eigen"
#include <array>
#include <concepts>

#include "Gedim_Macro.hpp"

#define DOFSMANAGER_MAX_DIMENSION 3

namespace Polydim
{
namespace PDETools
{
namespace DOFs
{
#if PYBIND == 1

#else
template <class mesh_connectivity_class>
concept is_mesh_connectivity_class_0D = requires(mesh_connectivity_class mesh) {
    { mesh.Cell0Ds_number() } -> std::same_as<unsigned int>;
    { mesh.Cell0D_marker(0) } -> std::same_as<unsigned int>;
};

template <class mesh_connectivity_class>
concept is_mesh_connectivity_class_1D = requires(mesh_connectivity_class mesh) {
    { mesh.Cell1Ds_number() } -> std::same_as<unsigned int>;
    { mesh.Cell1D_marker(0) } -> std::same_as<unsigned int>;
    { mesh.Cell1D_vertices(0) } -> std::same_as<std::array<unsigned int, 2>>;
};

template <class mesh_connectivity_class>
concept is_mesh_connectivity_class_2D = requires(mesh_connectivity_class mesh) {
    { mesh.Cell2Ds_number() } -> std::same_as<unsigned int>;
    { mesh.Cell2D_marker(0) } -> std::same_as<unsigned int>;
    { mesh.Cell2D_vertices(0) } -> std::same_as<std::vector<unsigned int>>;
    { mesh.Cell2D_edges(0) } -> std::same_as<std::vector<unsigned int>>;
};

template <class mesh_connectivity_class>
concept is_mesh_connectivity_class_3D = requires(mesh_connectivity_class mesh) {
    { mesh.Cell3Ds_number() } -> std::same_as<unsigned int>;
    { mesh.Cell3D_marker(0) } -> std::same_as<unsigned int>;
    { mesh.Cell3D_vertices(0) } -> std::same_as<std::vector<unsigned int>>;
    { mesh.Cell3D_edges(0) } -> std::same_as<std::vector<unsigned int>>;
    { mesh.Cell3D_faces(0) } -> std::same_as<std::vector<unsigned int>>;
};
#endif

class DOFsManager
{
  public:
    struct MeshDOFsInfo final
    {
        struct BoundaryInfo
        {
            enum struct BoundaryTypes
            {
                Unknwon = 0,
                Strong = 1,
                Weak = 2,
                None = 3
            };

            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes Type;
            unsigned int Marker;
        };

        std::array<std::vector<unsigned int>, DOFSMANAGER_MAX_DIMENSION + 1> CellsNumDOFs;
        std::array<std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, DOFSMANAGER_MAX_DIMENSION + 1> CellsBoundaryInfo;
    };

    using BoundaryTypes = typename MeshDOFsInfo::BoundaryInfo::BoundaryTypes;

    struct ConstantDOFsInfo final
    {
        std::array<unsigned int, DOFSMANAGER_MAX_DIMENSION + 1> NumDOFs;
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> BoundaryInfo;
    };

    struct DOFsData final
    {
        struct DOF final
        {
            enum struct Types
            {
                Unknwon = 0,
                Strong = 1,
                DOF = 2
            };

            Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types Type;
            unsigned int Global_Index;
        };

        struct GlobalCell_DOF
        {
            unsigned int Dimension;
            unsigned int CellIndex;
            unsigned int DOFIndex;
        };

        unsigned int NumberDOFs;
        unsigned int NumberInternalDOFs;
        unsigned int NumberBoundaryDOFs;
        unsigned int NumberStrongs;
        std::array<std::vector<std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF>>, DOFSMANAGER_MAX_DIMENSION + 1> CellsDOFs;
        std::array<std::vector<std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData::GlobalCell_DOF>>, DOFSMANAGER_MAX_DIMENSION + 1> CellsGlobalDOFs;
    };

    struct CellsDOFsIndicesData final
    {
        std::vector<std::vector<unsigned int>> Cells_DOFs_LocalIndex;
        std::vector<std::vector<unsigned int>> Cells_Strongs_LocalIndex;
        std::vector<std::vector<unsigned int>> Cells_DOFs_GlobalIndex;
        std::vector<std::vector<unsigned int>> Cells_Strongs_GlobalIndex;
    };

  private:
    void ConcatenateGlobalDOFs(const unsigned int local_cell_dimension,
                               const unsigned int local_cell_index,
                               const std::vector<typename DOFsData::DOF> &local_cell_DOFs,
                               std::vector<typename DOFsData::GlobalCell_DOF> &global_cell_DOFs,
                               unsigned int &globalDOF_counter) const;

    void CreateCellDOFs(const MeshDOFsInfo &meshDOFsInfo, DOFsData &dofs, const unsigned int dim) const;

    template <class mesh_connectivity_data_class>
    void Create_Constant_DOFsInfo_0D(const mesh_connectivity_data_class &mesh,
                                     const ConstantDOFsInfo &boundary_info,
                                     MeshDOFsInfo &mesh_dof_info) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_0D<mesh_connectivity_data_class>)
#endif
    {
        mesh_dof_info.CellsNumDOFs[0].resize(mesh.Cell0Ds_number(), boundary_info.NumDOFs[0]);
        mesh_dof_info.CellsBoundaryInfo[0].resize(mesh.Cell0Ds_number());

        for (unsigned int c = 0; c < mesh.Cell0Ds_number(); ++c)
        {
            mesh_dof_info.CellsBoundaryInfo[0][c] = boundary_info.BoundaryInfo.at(mesh.Cell0D_marker(c));
        }
    }

    template <class mesh_connectivity_data_class>
    void Create_Constant_DOFsInfo_1D(const mesh_connectivity_data_class &mesh,
                                     const ConstantDOFsInfo &boundary_info,
                                     MeshDOFsInfo &mesh_dof_info) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_1D<mesh_connectivity_data_class>)
#endif
    {
        mesh_dof_info.CellsNumDOFs[1].resize(mesh.Cell1Ds_number(), boundary_info.NumDOFs[1]);
        mesh_dof_info.CellsBoundaryInfo[1].resize(mesh.Cell1Ds_number());

        for (unsigned int c = 0; c < mesh.Cell1Ds_number(); ++c)
        {
            mesh_dof_info.CellsBoundaryInfo[1][c] = boundary_info.BoundaryInfo.at(mesh.Cell1D_marker(c));
        }
    }

    template <class mesh_connectivity_data_class>
    void Create_Constant_DOFsInfo_2D(const mesh_connectivity_data_class &mesh,
                                     const ConstantDOFsInfo &boundary_info,
                                     MeshDOFsInfo &mesh_dof_info) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_2D<mesh_connectivity_data_class>)
#endif
    {
        mesh_dof_info.CellsNumDOFs[2].resize(mesh.Cell2Ds_number(), boundary_info.NumDOFs[2]);
        mesh_dof_info.CellsBoundaryInfo[2].resize(mesh.Cell2Ds_number());

        for (unsigned int c = 0; c < mesh.Cell2Ds_number(); ++c)
        {
            mesh_dof_info.CellsBoundaryInfo[2][c] = boundary_info.BoundaryInfo.at(mesh.Cell2D_marker(c));
        }
    }

    template <class mesh_connectivity_data_class>
    void Create_Constant_DOFsInfo_3D(const mesh_connectivity_data_class &mesh,
                                     const ConstantDOFsInfo &boundary_info,
                                     MeshDOFsInfo &mesh_dof_info) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_3D<mesh_connectivity_data_class>)
#endif
    {
        mesh_dof_info.CellsNumDOFs[3].resize(mesh.Cell3Ds_number(), boundary_info.NumDOFs[3]);
        mesh_dof_info.CellsBoundaryInfo[3].resize(mesh.Cell3Ds_number());

        for (unsigned int c = 0; c < mesh.Cell3Ds_number(); ++c)
        {
            mesh_dof_info.CellsBoundaryInfo[3][c] = boundary_info.BoundaryInfo.at(mesh.Cell3D_marker(c));
        }
    }

#if PYBIND == 1
#else
    template <unsigned int dimension>
#endif
    void CreateCell0DDOFs(const MeshDOFsInfo &meshDOFsInfo, DOFsData &dofs) const
    {
#if PYBIND == 1
#else
        if (dimension < 0)
            return;
#endif

        CreateCellDOFs(meshDOFsInfo, dofs, 0);

        const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(0).size();

        auto &cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(0);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell0DIndex = 0; cell0DIndex < numCells; cell0DIndex++)
        {
            const auto &cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell0DIndex);

            const unsigned int cellNumGlobalDOFs = cell0D_DOFs.size();
            cellsGlobalDOFs[cell0DIndex].resize(cellNumGlobalDOFs);

            unsigned int globalDOF_counter = 0;

            ConcatenateGlobalDOFs(0, cell0DIndex, cell0D_DOFs, cellsGlobalDOFs[cell0DIndex], globalDOF_counter);
        }
    }

#if PYBIND == 1
    template <class mesh_connectivity_data_class>
#else
    template <unsigned int dimension, class mesh_connectivity_data_class>
#endif
    void CreateCell1DDOFs(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh, DOFsData &dofs) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_1D<mesh_connectivity_data_class>)
#endif
    {
#if PYBIND == 1
#else
        if (dimension < 1)
            return;
#endif

        CreateCellDOFs(meshDOFsInfo, dofs, 1);

        const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(1).size();

        auto &cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(1);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell1DIndex = 0; cell1DIndex < numCells; cell1DIndex++)
        {
            const auto cell1D_vertices = mesh.Cell1D_vertices(cell1DIndex);

            const unsigned int cell1D_origin_cell0DIndex = cell1D_vertices.at(0);
            const unsigned int cell1D_end_cell0DIndex = cell1D_vertices.at(1);

            const auto &origin_cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell1D_origin_cell0DIndex);
            const auto &end_cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell1D_end_cell0DIndex);
            const auto &cell1D_DOFs = dofs.CellsDOFs.at(1).at(cell1DIndex);

            const unsigned int cellNumGlobalDOFs = origin_cell0D_DOFs.size() + end_cell0D_DOFs.size() + cell1D_DOFs.size();
            cellsGlobalDOFs[cell1DIndex].resize(cellNumGlobalDOFs);

            unsigned int globalDOF_counter = 0;

            ConcatenateGlobalDOFs(0, cell1D_origin_cell0DIndex, origin_cell0D_DOFs, cellsGlobalDOFs[cell1DIndex], globalDOF_counter);
            ConcatenateGlobalDOFs(0, cell1D_end_cell0DIndex, end_cell0D_DOFs, cellsGlobalDOFs[cell1DIndex], globalDOF_counter);
            ConcatenateGlobalDOFs(1, cell1DIndex, cell1D_DOFs, cellsGlobalDOFs[cell1DIndex], globalDOF_counter);
        }
    }

#if PYBIND == 1
    template <class mesh_connectivity_data_class>
#else
    template <unsigned int dimension, class mesh_connectivity_data_class>
#endif
    void CreateCell2DDOFs(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh, DOFsData &dofs) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_2D<mesh_connectivity_data_class>)
#endif
    {
#if PYBIND == 1
#else
        if (dimension < 2)
            return;
#endif

        CreateCellDOFs(meshDOFsInfo, dofs, 2);

        const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(2).size();

        auto &cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(2);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell2DIndex = 0; cell2DIndex < numCells; cell2DIndex++)
        {
            unsigned int cellNumGlobalDOFs = dofs.CellsDOFs.at(2).at(cell2DIndex).size();

            const auto cell2D_vertices = mesh.Cell2D_vertices(cell2DIndex);

            for (unsigned int v = 0; v < cell2D_vertices.size(); v++)
            {
                const unsigned int vertex_cell0DIndex = cell2D_vertices.at(v);
                cellNumGlobalDOFs += dofs.CellsDOFs.at(0).at(vertex_cell0DIndex).size();
            }

            const auto cell2D_edges = mesh.Cell2D_edges(cell2DIndex);

            for (unsigned int e = 0; e < cell2D_edges.size(); e++)
            {
                const unsigned int edge_cell1DIndex = cell2D_edges.at(e);
                cellNumGlobalDOFs += dofs.CellsDOFs.at(1).at(edge_cell1DIndex).size();
            }

            cellsGlobalDOFs[cell2DIndex].resize(cellNumGlobalDOFs);

            unsigned int globalDOF_counter = 0;

            for (unsigned int v = 0; v < cell2D_vertices.size(); v++)
            {
                const unsigned int vertex_cell0DIndex = cell2D_vertices.at(v);
                ConcatenateGlobalDOFs(0, vertex_cell0DIndex, dofs.CellsDOFs.at(0).at(vertex_cell0DIndex), cellsGlobalDOFs[cell2DIndex], globalDOF_counter);
            }

            for (unsigned int e = 0; e < cell2D_edges.size(); e++)
            {
                const unsigned int edge_cell1DIndex = cell2D_edges.at(e);
                ConcatenateGlobalDOFs(1, edge_cell1DIndex, dofs.CellsDOFs.at(1).at(edge_cell1DIndex), cellsGlobalDOFs[cell2DIndex], globalDOF_counter);
            }

            ConcatenateGlobalDOFs(2, cell2DIndex, dofs.CellsDOFs.at(2).at(cell2DIndex), cellsGlobalDOFs[cell2DIndex], globalDOF_counter);
        }
    }

#if PYBIND == 1
    template <class mesh_connectivity_data_class>
#else
    template <unsigned int dimension, class mesh_connectivity_data_class>
#endif
    void CreateCell3DDOFs(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh, DOFsData &dofs) const
#if PYBIND == 1
#else
        requires(is_mesh_connectivity_class_3D<mesh_connectivity_data_class>)
#endif
    {
#if PYBIND == 1
#else
        if (dimension < 3)
            return;
#endif

        CreateCellDOFs(meshDOFsInfo, dofs, 3);

        const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(3).size();

        auto &cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(3);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell3DIndex = 0; cell3DIndex < numCells; cell3DIndex++)
        {
            unsigned int cellNumGlobalDOFs = dofs.CellsDOFs.at(3).at(cell3DIndex).size();
            const auto cell3D_vertices = mesh.Cell3D_vertices(cell3DIndex);
            const auto cell3D_edges = mesh.Cell3D_edges(cell3DIndex);
            const auto cell3D_faces = mesh.Cell3D_faces(cell3DIndex);

            for (unsigned int v = 0; v < cell3D_vertices.size(); v++)
            {
                const unsigned int vertex_cell0DIndex = cell3D_vertices.at(v);
                cellNumGlobalDOFs += dofs.CellsDOFs.at(0).at(vertex_cell0DIndex).size();
            }

            for (unsigned int e = 0; e < cell3D_edges.size(); e++)
            {
                const unsigned int edge_cell1DIndex = cell3D_edges.at(e);
                cellNumGlobalDOFs += dofs.CellsDOFs.at(1).at(edge_cell1DIndex).size();
            }

            for (unsigned int f = 0; f < cell3D_faces.size(); f++)
            {
                const unsigned int face_cell2DIndex = cell3D_faces.at(f);
                cellNumGlobalDOFs += dofs.CellsDOFs.at(2).at(face_cell2DIndex).size();
            }

            cellsGlobalDOFs[cell3DIndex].resize(cellNumGlobalDOFs);

            unsigned int globalDOF_counter = 0;

            for (unsigned int v = 0; v < cell3D_vertices.size(); v++)
            {
                const unsigned int vertex_cell0DIndex = cell3D_vertices.at(v);
                ConcatenateGlobalDOFs(0, vertex_cell0DIndex, dofs.CellsDOFs.at(0).at(vertex_cell0DIndex), cellsGlobalDOFs[cell3DIndex], globalDOF_counter);
            }

            for (unsigned int e = 0; e < cell3D_edges.size(); e++)
            {
                const unsigned int edge_cell1DIndex = cell3D_edges.at(e);
                ConcatenateGlobalDOFs(1, edge_cell1DIndex, dofs.CellsDOFs.at(1).at(edge_cell1DIndex), cellsGlobalDOFs[cell3DIndex], globalDOF_counter);
            }

            for (unsigned int f = 0; f < cell3D_faces.size(); f++)
            {
                const unsigned int face_cell2DIndex = cell3D_faces.at(f);
                ConcatenateGlobalDOFs(2, face_cell2DIndex, dofs.CellsDOFs.at(2).at(face_cell2DIndex), cellsGlobalDOFs[cell3DIndex], globalDOF_counter);
            }

            ConcatenateGlobalDOFs(3, cell3DIndex, dofs.CellsDOFs.at(3).at(cell3DIndex), cellsGlobalDOFs[cell3DIndex], globalDOF_counter);
        }
    }

  public:
    template <class mesh_connectivity_data_class>
    MeshDOFsInfo Create_Constant_DOFsInfo_0D(const mesh_connectivity_data_class &mesh, const ConstantDOFsInfo &boundary_info) const
    {
        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo meshDOFsInfo;

#if PYBIND == 1
        Create_Constant_DOFsInfo_0D(mesh, boundary_info, meshDOFsInfo);
#else
        Create_Constant_DOFsInfo_0D<0>(mesh, boundary_info, meshDOFsInfo);
#endif

        return meshDOFsInfo;
    }

    template <class mesh_connectivity_data_class>
    MeshDOFsInfo Create_Constant_DOFsInfo_1D(const mesh_connectivity_data_class &mesh, const ConstantDOFsInfo &boundary_info) const
    {
        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo meshDOFsInfo;

#if PYBIND == 1
        Create_Constant_DOFsInfo_0D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D(mesh, boundary_info, meshDOFsInfo);
#else
        Create_Constant_DOFsInfo_0D<1>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D<1>(mesh, boundary_info, meshDOFsInfo);
#endif

        return meshDOFsInfo;
    }

    template <class mesh_connectivity_data_class>
    MeshDOFsInfo Create_Constant_DOFsInfo_2D(const mesh_connectivity_data_class &mesh, const ConstantDOFsInfo &boundary_info) const
    {
        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo meshDOFsInfo;
#if PYBIND == 1
        Create_Constant_DOFsInfo_0D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_2D(mesh, boundary_info, meshDOFsInfo);
#else
        Create_Constant_DOFsInfo_0D<2>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D<2>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_2D<2>(mesh, boundary_info, meshDOFsInfo);
#endif
        return meshDOFsInfo;
    }

    template <class mesh_connectivity_data_class>
    MeshDOFsInfo Create_Constant_DOFsInfo_3D(const mesh_connectivity_data_class &mesh, const ConstantDOFsInfo &boundary_info) const
    {
        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo meshDOFsInfo;
#if PYBIND == 1
        Create_Constant_DOFsInfo_0D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_2D(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_3D(mesh, boundary_info, meshDOFsInfo);
#else
        Create_Constant_DOFsInfo_0D<3>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_1D<3>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_2D<3>(mesh, boundary_info, meshDOFsInfo);
        Create_Constant_DOFsInfo_3D<3>(mesh, boundary_info, meshDOFsInfo);
#endif

        return meshDOFsInfo;
    }

    DOFsData CreateDOFs_0D(const MeshDOFsInfo &meshDOFsInfo) const
    {
        DOFsData result;

        result.NumberDOFs = 0;
        result.NumberStrongs = 0;
        result.NumberBoundaryDOFs = 0;
        result.NumberInternalDOFs = 0;

#if PYBIND == 1
        CreateCell0DDOFs(meshDOFsInfo, result);
#else
        CreateCell0DDOFs<0>(meshDOFsInfo, result);
#endif

        return result;
    }

    template <class mesh_connectivity_data_class>
    DOFsData CreateDOFs_1D(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh) const
    {
        DOFsData result;

        result.NumberDOFs = 0;
        result.NumberStrongs = 0;
        result.NumberBoundaryDOFs = 0;
        result.NumberInternalDOFs = 0;

#if PYBIND == 1
        CreateCell0DDOFs(meshDOFsInfo, result);
        CreateCell1DDOFs(meshDOFsInfo, mesh, result);
#else
        CreateCell0DDOFs<1>(meshDOFsInfo, result);
        CreateCell1DDOFs<1>(meshDOFsInfo, mesh, result);
#endif

        return result;
    }

    template <class mesh_connectivity_data_class>
    DOFsData CreateDOFs_2D(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh) const
    {
        DOFsData result;

        result.NumberDOFs = 0;
        result.NumberStrongs = 0;
        result.NumberBoundaryDOFs = 0;
        result.NumberInternalDOFs = 0;

#if PYBIND == 1
        CreateCell0DDOFs(meshDOFsInfo, result);
        CreateCell1DDOFs(meshDOFsInfo, mesh, result);
        CreateCell2DDOFs(meshDOFsInfo, mesh, result);
#else
        CreateCell0DDOFs<2>(meshDOFsInfo, result);
        CreateCell1DDOFs<2>(meshDOFsInfo, mesh, result);
        CreateCell2DDOFs<2>(meshDOFsInfo, mesh, result);
#endif

        return result;
    }

    template <class mesh_connectivity_data_class>
    DOFsData CreateDOFs_3D(const MeshDOFsInfo &meshDOFsInfo, const mesh_connectivity_data_class &mesh) const
    {
        DOFsData result;

        result.NumberDOFs = 0;
        result.NumberStrongs = 0;
        result.NumberBoundaryDOFs = 0;
        result.NumberInternalDOFs = 0;
#if PYBIND == 1
        CreateCell0DDOFs(meshDOFsInfo, result);
        CreateCell1DDOFs(meshDOFsInfo, mesh, result);
        CreateCell2DDOFs(meshDOFsInfo, mesh, result);
        CreateCell3DDOFs(meshDOFsInfo, mesh, result);
#else
        CreateCell0DDOFs<3>(meshDOFsInfo, result);
        CreateCell1DDOFs<3>(meshDOFsInfo, mesh, result);
        CreateCell2DDOFs<3>(meshDOFsInfo, mesh, result);
        CreateCell3DDOFs<3>(meshDOFsInfo, mesh, result);
#endif

        return result;
    }

    CellsDOFsIndicesData ComputeCellsDOFsIndices(const DOFsData &dofs, const unsigned int dim) const;
};
} // namespace DOFs
} // namespace PDETools
} // namespace Polydim

#endif
