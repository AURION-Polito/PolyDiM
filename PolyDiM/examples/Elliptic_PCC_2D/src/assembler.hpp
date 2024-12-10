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
      void ComputeStrongTerm(const Gedim::GeometryUtilities& geometryUtilities,
                             const Gedim::MeshMatricesDAO& mesh,
                             const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                             const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                             const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                             const std::function<Eigen::VectorXd(const unsigned int,
                                                                 const Eigen::MatrixXd&)>& strong_boundary_condition,
                             Elliptic_PCC_2D_Problem_Data& assembler_data) const;

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
                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& diffusion_term,
                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& source_term,
                                            const std::function<Eigen::VectorXd(const unsigned int,
                                                                                const Eigen::MatrixXd&)>& strong_boundary_condition) const;

      VEM_Performance_Result ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                   const Gedim::MeshMatricesDAO& mesh,
                                                   const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                   const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data) const;

      PostProcess_Data PostProcessSolution(const Gedim::GeometryUtilities& geometryUtilities,
                                           const Gedim::MeshMatricesDAO& mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                           const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                           const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                           const Elliptic_PCC_2D_Problem_Data& assembler_data,
                                           const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& exact_solution,
                                           const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& exact_derivative_solution) const;
  };

}

#endif
