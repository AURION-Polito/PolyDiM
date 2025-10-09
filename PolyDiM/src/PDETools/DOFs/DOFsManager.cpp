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

#include "DOFsManager.hpp"

namespace Polydim
{
  namespace PDETools
  {
    namespace DOFs
    {
      //***************************************************************************
      void DOFsManager::ConcatenateGlobalDOFs(const unsigned int local_cell_dimension, const unsigned int local_cell_index, const std::vector<DOFsData::DOF>& local_cell_DOFs, std::vector<DOFsData::GlobalCell_DOF>& global_cell_DOFs, unsigned int& globalDOF_counter) const
      {
        for (unsigned int d = 0; d < local_cell_DOFs.size(); d++)
        {
          global_cell_DOFs[globalDOF_counter] = {local_cell_dimension, local_cell_index, d};

          globalDOF_counter++;
        }
      }
      //***************************************************************************
      void DOFsManager::CreateCellDOFs(const MeshDOFsInfo& meshDOFsInfo, DOFsData& dofs, const unsigned int dim) const
      {
        const auto &cells_num_dofs = meshDOFsInfo.CellsNumDOFs.at(dim);
        const auto &cells_boundary_info = meshDOFsInfo.CellsBoundaryInfo.at(dim);
        const unsigned int numCells = cells_num_dofs.size();
        auto &cellsDOFs = dofs.CellsDOFs.at(dim);
        cellsDOFs.resize(numCells);

        for (unsigned int c = 0; c < numCells; c++)
        {
          const unsigned int numCellDofs = cells_num_dofs.at(c);

          const auto &cell_boundary_info = cells_boundary_info.at(c);
          const BoundaryTypes &cellBoundaryType = cell_boundary_info.Type;

          cellsDOFs.at(c).resize(numCellDofs);

          switch (cellBoundaryType)
          {
            case BoundaryTypes::None: {
              for (unsigned int d = 0; d < numCellDofs; d++)
              {
                auto &dof = cellsDOFs.at(c).at(d);
                dof.Type = DOFsData::DOF::Types::DOF;
                dof.Global_Index = dofs.NumberDOFs + d;
              }

              dofs.NumberDOFs += numCellDofs;
              dofs.NumberInternalDOFs += numCellDofs;
            }
              break;
            case BoundaryTypes::Strong: {
              for (unsigned int d = 0; d < numCellDofs; d++)
              {
                auto &dof = cellsDOFs.at(c).at(d);
                dof.Type = DOFsData::DOF::Types::Strong;
                dof.Global_Index = dofs.NumberStrongs + d;
              }

              dofs.NumberStrongs += numCellDofs;
            }
              break;
            case BoundaryTypes::Weak: {
              for (unsigned int d = 0; d < numCellDofs; d++)
              {
                auto &dof = cellsDOFs.at(c).at(d);
                dof.Type = DOFsData::DOF::Types::DOF;
                dof.Global_Index = dofs.NumberDOFs + d;
              }

              dofs.NumberDOFs += numCellDofs;
              dofs.NumberBoundaryDOFs += numCellDofs;
            }
              break;
            default:
              throw std::runtime_error("Unknown BoundaryTypes");
              break;
          }
        }
      }
      //***************************************************************************
      DOFsManager::CellsDOFsIndicesData DOFsManager::ComputeCellsDOFsIndices(const DOFsData& dofs, const unsigned int dimension) const
      {
        CellsDOFsIndicesData result;

        const unsigned int num_cells = dofs.CellsGlobalDOFs[dimension].size();

        result.Cells_DOFs_LocalIndex.resize(num_cells);
        result.Cells_DOFs_GlobalIndex.resize(num_cells);
        result.Cells_Strongs_LocalIndex.resize(num_cells);
        result.Cells_Strongs_GlobalIndex.resize(num_cells);

        for (unsigned int c = 0; c < num_cells; ++c)
        {
          const auto& global_dofs = dofs.CellsGlobalDOFs[0][c];

          std::list<unsigned int> dofs_local_index;
          std::list<unsigned int> dofs_global_index;
          std::list<unsigned int> strongs_local_index;
          std::list<unsigned int> strongs_global_index;

          for (unsigned int g_d = 0; g_d < global_dofs.size(); ++g_d)
          {
            const auto& global_dof = global_dofs.at(g_d);
            const auto &local_dof = dofs.CellsDOFs.at(global_dof.Dimension).at(global_dof.CellIndex).at(global_dof.DOFIndex);

            switch (local_dof.Type)
            {
              case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
              {
                strongs_local_index.push_back(g_d);
                strongs_global_index.push_back(local_dof.Global_Index);
              }
                break;
              case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
              {
                dofs_local_index.push_back(g_d);
                dofs_global_index.push_back(local_dof.Global_Index);
              }
                break;
              default:
                throw std::runtime_error("Unknown DOF Type");
            }
          }

          result.Cells_DOFs_LocalIndex[c] = std::vector<unsigned int>(dofs_local_index.begin(),
                                                                      dofs_local_index.end());
          result.Cells_DOFs_GlobalIndex[c] = std::vector<unsigned int>(dofs_global_index.begin(),
                                                                       dofs_global_index.end());
          result.Cells_Strongs_LocalIndex[c] = std::vector<unsigned int>(strongs_local_index.begin(),
                                                                         strongs_local_index.end());
          result.Cells_Strongs_GlobalIndex[c] = std::vector<unsigned int>(strongs_global_index.begin(),
                                                                          strongs_global_index.end());
        }

        return result;
      }
      //***************************************************************************
    } // namespace DOFs
  } // namespace PDETools
} // namespace Polydim
