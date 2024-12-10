#ifndef __PDETOOLS_DOFS_DOFsManager_HPP
#define __PDETOOLS_DOFS_DOFsManager_HPP

#include "Eigen/Eigen"
#include <concepts>
#include <array>

namespace Polydim
{
  namespace PDETools
  {
    namespace DOFs
    {
      template <class mesh_connectivity_class>
      concept is_mesh_connectivity_class_1D = requires(mesh_connectivity_class mesh) {
      { mesh.Cell1D_vertices(0) } -> std::same_as<std::array<unsigned int, 2>>; };

      template <class mesh_connectivity_class>
      concept is_mesh_connectivity_class_2D = requires(mesh_connectivity_class mesh) {
      { mesh.Cell2D_vertices(0) } -> std::same_as<std::vector<unsigned int>>;
      { mesh.Cell2D_edges(0) } -> std::same_as<std::vector<unsigned int>>; };

      template <class mesh_connectivity_class>
      concept is_mesh_connectivity_class_3D = requires(mesh_connectivity_class mesh) {
      { mesh.Cell3D_vertices(0) } -> std::same_as<std::vector<unsigned int>>;
      { mesh.Cell3D_edges(0) } -> std::same_as<std::vector<unsigned int>>;
      { mesh.Cell3D_faces(0) } -> std::same_as<std::vector<unsigned int>>; };

      template <unsigned int dimension>
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

                  BoundaryTypes Type;
                  unsigned int Marker;
              };

              std::array<std::vector<unsigned int>, dimension + 1> CellsNumDOFs;
              std::array<std::vector<BoundaryInfo>, dimension + 1> CellsBoundaryInfo;
          };

          using BoundaryTypes = typename MeshDOFsInfo::BoundaryInfo::BoundaryTypes;

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

                  Types Type;
                  unsigned int Global_Index;
              };

              struct GlobalCell_DOF
              {
                  unsigned int Dimension;
                  unsigned int CellIndex;
                  unsigned int DOFIndex;
              };

              unsigned int NumberDOFs;
              unsigned int NumberStrongs;
              std::array<std::vector<std::vector<DOF>>, dimension + 1> CellsDOFs;
              std::array<std::vector<std::vector<GlobalCell_DOF>>, dimension + 1> CellsGlobalDOFs;
          };

        private:
          void ConcatenateGlobalDOFs(const unsigned int local_cell_dimension,
                                     const unsigned int local_cell_index,
                                     const std::vector<typename DOFsData::DOF>& local_cell_DOFs,
                                     std::vector<typename DOFsData::GlobalCell_DOF>& global_cell_DOFs,
                                     unsigned int& globalDOF_counter) const
          {
            for (unsigned int d = 0; d < local_cell_DOFs.size(); d++)
            {
              global_cell_DOFs[globalDOF_counter] =
              {
                local_cell_dimension,
                local_cell_index,
                d
              };

              globalDOF_counter++;
            }
          }

          void CreateCellDOFs(const MeshDOFsInfo& meshDOFsInfo,
                              DOFsData& dofs,
                              const unsigned int dim) const
          {
            const auto& cells_num_dofs = meshDOFsInfo.CellsNumDOFs.at(dim);
            const auto& cells_boundary_info = meshDOFsInfo.CellsBoundaryInfo.at(dim);
            const unsigned int numCells = cells_num_dofs.size();
            auto& cellsDOFs = dofs.CellsDOFs.at(dim);
            cellsDOFs.resize(numCells);

            for (unsigned int c = 0; c < numCells; c++)
            {
              const unsigned int numCellDofs = cells_num_dofs.at(c);

              const auto& cell_boundary_info = cells_boundary_info.at(c);
              const BoundaryTypes& cellBoundaryType = cell_boundary_info.Type;

              cellsDOFs.at(c).resize(numCellDofs);

              switch (cellBoundaryType)
              {
                case BoundaryTypes::None:
                {
                  for (unsigned int d = 0; d < numCellDofs; d++)
                  {
                    auto& dof = cellsDOFs.at(c).at(d);
                    dof.Type = DOFsData::DOF::Types::DOF;
                    dof.Global_Index = dofs.NumberDOFs + d;
                  }

                  dofs.NumberDOFs += numCellDofs;
                }
                  break;
                case BoundaryTypes::Strong:
                {
                  for (unsigned int d = 0; d < numCellDofs; d++)
                  {
                    auto& dof = cellsDOFs.at(c).at(d);
                    dof.Type = DOFsData::DOF::Types::Strong;
                    dof.Global_Index = dofs.NumberStrongs + d;
                  }

                  dofs.NumberStrongs += numCellDofs;
                }
                  break;
                case BoundaryTypes::Weak:
                {
                  for (unsigned int d = 0; d < numCellDofs; d++)
                  {
                    auto& dof = cellsDOFs.at(c).at(d);
                    dof.Type = DOFsData::DOF::Types::DOF;
                    dof.Global_Index = dofs.NumberDOFs + d;
                  }

                  dofs.NumberDOFs += numCellDofs;
                }
                  break;
                default:
                  throw std::runtime_error("Unknown BoundaryTypes");
                  break;
              }
            }
          }

          void CreateCell0DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                DOFsData& dofs) const
          {
            if (dimension < 0)
              return;

            CreateCellDOFs(meshDOFsInfo,
                           dofs,
                           0);

            const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(0).size();

            auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(0);
            cellsGlobalDOFs.resize(numCells);

            for (unsigned int cell0DIndex = 0; cell0DIndex < numCells; cell0DIndex++)
            {
              const auto& cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell0DIndex);

              const unsigned int cellNumGlobalDOFs = cell0D_DOFs.size();
              cellsGlobalDOFs[cell0DIndex].resize(cellNumGlobalDOFs);

              unsigned int globalDOF_counter = 0;

              ConcatenateGlobalDOFs(0,
                                    cell0DIndex,
                                    cell0D_DOFs,
                                    cellsGlobalDOFs[cell0DIndex],
                                    globalDOF_counter);
            }
          }

          template <class mesh_connectivity_data_class>
          void CreateCell1DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                const mesh_connectivity_data_class& mesh,
                                DOFsData& dofs) const requires(is_mesh_connectivity_class_1D<mesh_connectivity_data_class>)
          {
            if (dimension < 1)
              return;

            CreateCellDOFs(meshDOFsInfo,
                           dofs,
                           1);

            const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(1).size();

            auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(1);
            cellsGlobalDOFs.resize(numCells);

            for (unsigned int cell1DIndex = 0; cell1DIndex < numCells; cell1DIndex++)
            {
              const auto cell1D_vertices = mesh.Cell1D_vertices(cell1DIndex);

              const unsigned int cell1D_origin_cell0DIndex = cell1D_vertices.at(0);
              const unsigned int cell1D_end_cell0DIndex = cell1D_vertices.at(0);

              const auto& origin_cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell1D_origin_cell0DIndex);
              const auto& end_cell0D_DOFs = dofs.CellsDOFs.at(0).at(cell1D_end_cell0DIndex);
              const auto& cell1D_DOFs = dofs.CellsDOFs.at(1).at(cell1DIndex);

              const unsigned int cellNumGlobalDOFs = origin_cell0D_DOFs.size() +
                                                     end_cell0D_DOFs.size() +
                                                     cell1D_DOFs.size();
              cellsGlobalDOFs[cell1DIndex].resize(cellNumGlobalDOFs);

              unsigned int globalDOF_counter = 0;

              ConcatenateGlobalDOFs(0,
                                    cell1D_origin_cell0DIndex,
                                    origin_cell0D_DOFs,
                                    cellsGlobalDOFs[cell1DIndex],
                                    globalDOF_counter);
              ConcatenateGlobalDOFs(0,
                                    cell1D_end_cell0DIndex,
                                    end_cell0D_DOFs,
                                    cellsGlobalDOFs[cell1DIndex],
                                    globalDOF_counter);
              ConcatenateGlobalDOFs(1,
                                    cell1DIndex,
                                    cell1D_DOFs,
                                    cellsGlobalDOFs[cell1DIndex],
                                    globalDOF_counter);
            }
          }

          template <class mesh_connectivity_data_class>
          void CreateCell2DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                const mesh_connectivity_data_class& mesh,
                                DOFsData& dofs) const requires(is_mesh_connectivity_class_2D<mesh_connectivity_data_class>)
          {
            if (dimension < 2)
              return;

            CreateCellDOFs(meshDOFsInfo,
                           dofs,
                           2);

            const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(2).size();

            auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(2);
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
                ConcatenateGlobalDOFs(0,
                                      vertex_cell0DIndex,
                                      dofs.CellsDOFs.at(0).at(vertex_cell0DIndex),
                                      cellsGlobalDOFs[cell2DIndex],
                                      globalDOF_counter);
              }

              for (unsigned int e = 0; e < cell2D_edges.size(); e++)
              {
                const unsigned int edge_cell1DIndex = cell2D_edges.at(e);
                ConcatenateGlobalDOFs(1,
                                      edge_cell1DIndex,
                                      dofs.CellsDOFs.at(1).at(edge_cell1DIndex),
                                      cellsGlobalDOFs[cell2DIndex],
                                      globalDOF_counter);
              }

              ConcatenateGlobalDOFs(2,
                                    cell2DIndex,
                                    dofs.CellsDOFs.at(2).at(cell2DIndex),
                                    cellsGlobalDOFs[cell2DIndex],
                                    globalDOF_counter);
            }
          }

          template <class mesh_connectivity_data_class>
          void CreateCell3DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                const mesh_connectivity_data_class& mesh,
                                DOFsData& dofs) const requires(is_mesh_connectivity_class_3D<mesh_connectivity_data_class>)
          {
            if (dimension < 3)
              return;

            CreateCellDOFs(meshDOFsInfo,
                           dofs,
                           3);

            const unsigned int numCells = meshDOFsInfo.CellsNumDOFs.at(3).size();

            auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(3);
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
                ConcatenateGlobalDOFs(0,
                                      vertex_cell0DIndex,
                                      dofs.CellsDOFs.at(0).at(vertex_cell0DIndex),
                                      cellsGlobalDOFs[cell3DIndex],
                                      globalDOF_counter);
              }

              for (unsigned int e = 0; e < cell3D_edges.size(); e++)
              {
                const unsigned int edge_cell1DIndex = cell3D_edges.at(e);
                ConcatenateGlobalDOFs(1,
                                      edge_cell1DIndex,
                                      dofs.CellsDOFs.at(1).at(edge_cell1DIndex),
                                      cellsGlobalDOFs[cell3DIndex],
                                      globalDOF_counter);
              }

              for (unsigned int f = 0; f < cell3D_faces.size(); f++)
              {
                const unsigned int face_cell2DIndex = cell3D_faces.at(f);
                ConcatenateGlobalDOFs(2,
                                      face_cell2DIndex,
                                      dofs.CellsDOFs.at(2).at(face_cell2DIndex),
                                      cellsGlobalDOFs[cell3DIndex],
                                      globalDOF_counter);
              }

              ConcatenateGlobalDOFs(3,
                                    cell3DIndex,
                                    dofs.CellsDOFs.at(3).at(cell3DIndex),
                                    cellsGlobalDOFs[cell3DIndex],
                                    globalDOF_counter);
            }
          }

        public:
          DOFsData CreateDOFs(const MeshDOFsInfo& meshDOFsInfo) const requires(dimension == 0)
          {
            DOFsData result;

            result.NumberDOFs = 0;
            result.NumberStrongs = 0;

            CreateCell0DDOFs(meshDOFsInfo,
                             result);

            return result;
          }

          template <class mesh_connectivity_data_class>
          DOFsData CreateDOFs(const MeshDOFsInfo& meshDOFsInfo,
                              const mesh_connectivity_data_class& mesh) const requires(dimension == 1)
          {
            DOFsData result;

            result.NumberDOFs = 0;
            result.NumberStrongs = 0;

            CreateCell0DDOFs(meshDOFsInfo,
                             result);

            CreateCell1DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            return result;
          }

          template <class mesh_connectivity_data_class>
          DOFsData CreateDOFs(const MeshDOFsInfo& meshDOFsInfo,
                              const mesh_connectivity_data_class& mesh) const requires(dimension == 2)
          {
            DOFsData result;

            result.NumberDOFs = 0;
            result.NumberStrongs = 0;

            CreateCell0DDOFs(meshDOFsInfo,
                             result);

            CreateCell1DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            CreateCell2DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            return result;
          }

          template <class mesh_connectivity_data_class>
          DOFsData CreateDOFs(const MeshDOFsInfo& meshDOFsInfo,
                              const mesh_connectivity_data_class& mesh) const requires(dimension == 3)
          {
            DOFsData result;

            result.NumberDOFs = 0;
            result.NumberStrongs = 0;

            CreateCell0DDOFs(meshDOFsInfo,
                             result);

            CreateCell1DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            CreateCell2DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            CreateCell3DDOFs(meshDOFsInfo,
                             mesh,
                             result);

            return result;
          }
      };
    }
  }
}

#endif
