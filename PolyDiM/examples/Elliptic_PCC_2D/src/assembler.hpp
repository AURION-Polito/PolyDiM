#ifndef __assembler_H
#define __assembler_H

#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"
#include "Quadrature_Gauss1D.hpp"

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "VEM_PCC_Utilities.hpp"
#include "VEM_PCC_2D_Creator.hpp"
#include "DOFsManager.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"


namespace Polydim
{
  namespace examples
  {
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

          struct PostProcess_Data final
          {
              Eigen::VectorXd cell0Ds_numeric;
              Eigen::VectorXd cell0Ds_exact;

              Eigen::VectorXd cell2Ds_error_L2;
              Eigen::VectorXd cell2Ds_norm_L2;
              double error_L2;
              double norm_L2;
              Eigen::VectorXd cell2Ds_error_H1;
              Eigen::VectorXd cell2Ds_norm_H1;
              double error_H1;
              double norm_H1;

              double mesh_size;

              double residual_norm;
          };

        private:
          void ComputeStrongTerm(const Gedim::MeshMatricesDAO& mesh,
                                 const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                 const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                 const Polydim::examples::Elliptic_PCC_2D::test::I_Test& test,
                                 Elliptic_PCC_2D_Problem_Data& assembler_data) const;

          void ComputeWeakTerm(const unsigned int cell2DIndex,
                               const Gedim::MeshMatricesDAO& mesh,
                               const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry& polygon,
                               const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                               const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                               const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                               const Polydim::VEM::PCC::I_VEM_PCC_2D_LocalSpace& vem_local_space,
                               const Polydim::examples::Elliptic_PCC_2D::test::I_Test& test,
                               Elliptic_PCC_2D_Problem_Data& assembler_data) const;

        public:
          Elliptic_PCC_2D_Problem_Data Assemble(const Polydim::examples::Elliptic_PCC_2D::Program_configuration& config,
                                                const Gedim::MeshMatricesDAO& mesh,
                                                const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                                const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                const Polydim::examples::Elliptic_PCC_2D::test::I_Test& test) const;

          VEM_Performance_Result ComputeVemPerformance(const Polydim::examples::Elliptic_PCC_2D::Program_configuration& config,
                                                       const Gedim::MeshMatricesDAO& mesh,
                                                       const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                       const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data) const;

          PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_PCC_2D::Program_configuration& config,
                                               const Gedim::MeshMatricesDAO& mesh,
                                               const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                               const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                               const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                               const Elliptic_PCC_2D_Problem_Data& assembler_data,
                                               const Polydim::examples::Elliptic_PCC_2D::test::I_Test& test) const;

      };
    }
  }
}

#endif
