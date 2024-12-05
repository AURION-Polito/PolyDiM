#include "DOFsManager.hpp"
#include <iostream>

namespace Polydim
{
  namespace PDETools
  {
    namespace DOFs
    {
      // ***************************************************************************
      template class DOFsManager<0>;
      template class DOFsManager<1>;
      template class DOFsManager<2>;
      template class DOFsManager<3>;
      // ***************************************************************************
      template<unsigned int dimension>
      typename DOFsManager<dimension>::DOFsData DOFsManager<dimension>::CreateDOFs(const MeshDOFsInfo& meshDOFsInfo) const
      {
        DOFsData result;

        result.NumberDOFs = 0;
        result.NumberStrongs = 0;

        CreateCell0DDOFs(meshDOFsInfo,
                         result);

        CreateCell1DDOFs(meshDOFsInfo,
                         result);

        CreateCell2DDOFs(meshDOFsInfo,
                         result);

        CreateCell3DDOFs(meshDOFsInfo,
                         result);

        return result;
      }
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::ConcatenateGlobalDOFs(const unsigned int local_cell_dimension,
                                                         const unsigned int local_cell_index,
                                                         const std::vector<typename DOFsManager<dimension>::DOFsData::DOF>& local_cell_DOFs,
                                                         std::vector<typename DOFsManager<dimension>::DOFsData::GlobalCell_DOF>& global_cell_DOFs,
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
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::CreateCellDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                                  DOFsData& dofs,
                                                  const unsigned int d) const
      {
        const auto& cells_num_dofs = meshDOFsInfo.CellsNumDOFs.at(d);
        const auto& cells_boundary_info = meshDOFsInfo.CellsBoundaryInfo.at(d);
        const unsigned int numCells = cells_num_dofs.size();
        auto& cellsDOFs = dofs.CellsDOFs.at(d);
        cellsDOFs.resize(numCells);

        for (unsigned int c = 0; c < numCells; c++)
        {
          const unsigned int numCellDofs = cells_num_dofs.at(c);

          const auto& cell_boundary_info = cells_boundary_info.at(c);
          const unsigned int cellMarker = cell_boundary_info.Marker;
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
                dof.Boundary = { cellBoundaryType, cellMarker };
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
                dof.Boundary = { cellBoundaryType, cellMarker };
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
                dof.Boundary = { cellBoundaryType, cellMarker };
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
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::CreateCell0DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                                    DOFsData& dofs) const
      {
        if (dimension < 0)
          return;

        CreateCellDOFs(meshDOFsInfo,
                       dofs,
                       0);

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
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::CreateCell1DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                                    DOFsData& dofs) const
      {
        if (dimension < 1)
          return;

        CreateCellDOFs(meshDOFsInfo,
                       dofs,
                       1);


        auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(1);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell1DIndex = 0; cell1DIndex < numCells; cell1DIndex++)
        {
          const unsigned int cell1D_origin_cell0DIndex = mesh.Cell1Ds(0, cell1DIndex);
          const unsigned int cell1D_end_cell0DIndex = mesh.Cell1Ds(1, cell1DIndex);

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
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::CreateCell2DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                                    DOFsData& dofs) const
      {
        if (dimension < 2)
          return;

        CreateCellDOFs(meshDOFsInfo,
                       dofs,
                       2);

        auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(2);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell2DIndex = 0; cell2DIndex < numCells; cell2DIndex++)
        {
          const auto& cell2D = mesh.Cell2Ds.at(cell2DIndex);

          unsigned int cellNumGlobalDOFs = dofs.CellsDOFs.at(2).at(cell2DIndex).size();

          for (unsigned int v = 0; v < cell2D.cols(); v++)
          {
            const unsigned int vertex_cell0DIndex = cell2D(0, v);
            cellNumGlobalDOFs += dofs.CellsDOFs.at(0).at(vertex_cell0DIndex).size();
          }

          for (unsigned int e = 0; e < cell2D.cols(); e++)
          {
            const unsigned int edge_cell1DIndex = cell2D(1, e);
            cellNumGlobalDOFs += dofs.CellsDOFs.at(1).at(edge_cell1DIndex).size();
          }

          cellsGlobalDOFs[cell2DIndex].resize(cellNumGlobalDOFs);

          unsigned int globalDOF_counter = 0;

          for (unsigned int v = 0; v < cell2D.cols(); v++)
          {
            const unsigned int vertex_cell0DIndex = cell2D(0, v);
            ConcatenateGlobalDOFs(0,
                                  vertex_cell0DIndex,
                                  dofs.CellsDOFs.at(0).at(vertex_cell0DIndex),
                                  cellsGlobalDOFs[cell2DIndex],
                                  globalDOF_counter);
          }

          for (unsigned int e = 0; e < cell2D.cols(); e++)
          {
            const unsigned int edge_cell1DIndex = cell2D(1, e);
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
      // ***************************************************************************
      template<unsigned int dimension>
      void DOFsManager<dimension>::CreateCell3DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                                    DOFsData& dofs) const
      {
        if (dimension < 3)
          return;

        CreateCellDOFs(meshDOFsInfo,
                       dofs,
                       3);

        auto& cellsGlobalDOFs = dofs.CellsGlobalDOFs.at(3);
        cellsGlobalDOFs.resize(numCells);

        for (unsigned int cell3DIndex = 0; cell3DIndex < numCells; cell3DIndex++)
        {
          unsigned int cellNumGlobalDOFs = dofs.CellsDOFs.at(3).at(cell3DIndex).size();
          const auto cell3D_vertices = mesh.Cell3DsVertices[cell3DIndex];
          const auto cell3D_edges = mesh.Cell3DsEdges[cell3DIndex];
          const auto cell3D_faces = mesh.Cell3DsFaces[cell3DIndex];

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
      // ***************************************************************************
    }
  }
}
