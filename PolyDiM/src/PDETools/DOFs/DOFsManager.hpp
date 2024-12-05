#ifndef __PDETOOLS_DOFS_DOFsManager_HPP
#define __PDETOOLS_DOFS_DOFsManager_HPP

#include "Eigen/Eigen"

namespace Polydim
{
  namespace PDETools
  {
    namespace DOFs
    {
      struct DOFsManager_Mesh_Connectivity_Data final
      {
          Eigen::MatrixXi& Cell1Ds;
          std::vector<Eigen::MatrixXi>& Cell2Ds;
          std::vector<std::vector<unsigned int>>& Cell3DsVertices;
          std::vector<std::vector<unsigned int>>& Cell3DsEdges;
          std::vector<std::vector<unsigned int>>& Cell3DsFaces;
      };

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
                  typename MeshDOFsInfo::BoundaryInfo Boundary;
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
                                     const std::vector<typename DOFsManager<dimension>::DOFsData::DOF>& local_cell_DOFs,
                                     std::vector<typename DOFsManager<dimension>::DOFsData::GlobalCell_DOF>& global_cell_DOFs,
                                     unsigned int& globalDOF_counter) const;

          void CreateCellDOFs(const MeshDOFsInfo& meshDOFsInfo,
                              typename DOFsManager<dimension>::DOFsData& dofs,
                              const unsigned int d) const;

          void CreateCell0DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                typename DOFsManager<dimension>::DOFsData& dofs) const;
          void CreateCell1DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                typename DOFsManager<dimension>::DOFsData& dofs) const;
          void CreateCell2DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                typename DOFsManager<dimension>::DOFsData& dofs) const;
          void CreateCell3DDOFs(const MeshDOFsInfo& meshDOFsInfo,
                                typename DOFsManager<dimension>::DOFsData& dofs) const;

        public:
          typename DOFsManager<dimension>::DOFsData CreateDOFs(const MeshDOFsInfo& meshDOFsInfo) const;
      };
    }
  }
}

#endif
