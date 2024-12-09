#ifndef __assembler_H
#define __assembler_H

#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "DOFsManager.hpp"



namespace Elliptic_PCC_2D
{
  class Assembler final
  {
    public:

      struct Elliptic_PCC_2D_Problem_Data final
      {
          Gedim::Eigen_SparseArray<> globalMatrixA;
          Gedim::Eigen_SparseArray<> dirichletMatrixA;
          Gedim::Eigen_Array<> rightHandSide;
          Gedim::Eigen_Array<> solution;
          Gedim::Eigen_Array<> solutionDirichlet;
      };

      struct VEM_Performance_Result final
      {
          struct Cell2D_Performance final
          {
              unsigned int NumBoundaryQuadraturePoints = 0;
              unsigned int NumInternalQuadraturePoints = 0;
              Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
          };

          std::vector<Cell2D_Performance> Cell2DsPerformance;
      };

    private:
//      void ComputeStrongTerm(const Gedim::GeometryUtilities geometryUtilities,
//                             const Gedim::IMeshDAO& mesh,
//                             const Gedim::IDOFManagement& dofManager,
//                             const Gedim::IStrongBoundaryCondition& strongBoundaryCondition,
//                             Gedim::IArray& solutionDirichlet) const;

//      void ComputeWeakTerm(const Gedim::IMeshDAO& mesh,
//                           const unsigned int& cell2DIndex,
//                           const Eigen::VectorXd& cell2DEdgeLengths,
//                           const Eigen::MatrixXd& cell2DEdgeTangents,
//                           const Eigen::MatrixXd& cell2DEdgeNormals,
//                           const Eigen::MatrixXd& cell2DVertices,
//                           const std::vector<bool>& cell2DEdgeDirections,
//                           const Gedim::IDOFManagement& dofManager,
//                           const Gedim::VEM_IValues_PCC_2D& vemValues,
//                           const Gedim::VEM_ValuesData& vemLocalSpace,
//                           const Gedim::VEM_IQuadrature2D& vemQuadrature,
//                           const Gedim::IWeakBoundaryCondition& weakBoundaryCondition,
//                           Gedim::IArray& rightHandSide) const;
    public:
      Elliptic_PCC_2D_Problem_Data Assemble(const Gedim::GeometryUtilities& geometryUtilities,
                                            const Gedim::MeshMatricesDAO& mesh,
                                            const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                            const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                            const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& diffusionTerm,
                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& sourceTerm) const;

      VEM_Performance_Result ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                   const Gedim::IMeshDAO& mesh,
                                                   const std::vector<Eigen::MatrixXd>& meshCell2DsVertices,
                                                   const std::vector<std::vector<bool>>& meshCell2DsEdgeDirections,
                                                   const std::vector<std::vector<Eigen::Matrix3d> >& meshCell2DsTriangulations,
                                                   const std::vector<double>& meshCell2DsAreas,
                                                   const std::vector<Eigen::Vector3d>& meshCell2DsCentroids,
                                                   const std::vector<double>& meshCell2DsDiameters,
                                                   const std::vector<Eigen::VectorXd>& meshCell2DsEdgeLengths,
                                                   const std::vector<Eigen::MatrixXd>& meshCell2DsEdgeTangents,
                                                   const std::vector<Eigen::MatrixXd>& meshCell2DsEdgeNormals) const;
  };

}

#endif
