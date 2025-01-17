#ifndef __assembler_H
#define __assembler_H

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

#include "local_space.hpp"

namespace Polydim
{
namespace examples
{
namespace Elastic_PCC_2D
{
class Assembler final
{
  public:
    struct Elastic_PCC_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionDirichlet;
    };

    struct Performance_Data final
    {
        std::vector<local_space::Performance_Data> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        std::array<Eigen::VectorXd, 3> cell0Ds_numeric_displacement;
        std::array<Eigen::VectorXd, 3> cell0Ds_exact_displacement;

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
    void ComputeStrongTerm(const unsigned int cell2D_index,
                           const Gedim::MeshMatricesDAO &mesh,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                           const local_space::ReferenceElement_Data &reference_element_data,
                           const local_space::LocalSpace_Data &local_space_data,
                           const Polydim::examples::Elastic_PCC_2D::test::I_Test &test,
                           Elastic_PCC_2D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                         const local_space::ReferenceElement_Data &reference_element_data,
                         const local_space::LocalSpace_Data &local_space_data,
                         const Polydim::examples::Elastic_PCC_2D::test::I_Test &test,
                         Elastic_PCC_2D_Problem_Data &assembler_data) const;

    Eigen::MatrixXd ComputeSUPGMatrix(const std::array<Eigen::VectorXd, 3> &advection_term_values,
                                      const Eigen::VectorXd &diffusion_term_values,
                                      const Eigen::MatrixXd &basis_functions_values,
                                      const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                      const Eigen::VectorXd &quadrature_weights) const;

    Eigen::MatrixXd ComputeSUPGForcingTerm(const std::array<Eigen::VectorXd, 3> &advection_term_values,
                                           const Eigen::VectorXd &forcing_term_values,
                                           const std::vector<Eigen::MatrixXd> &basis_functions_derivative_values,
                                           const Eigen::VectorXd &quadrature_weights) const;

  public:
    Elastic_PCC_2D_Problem_Data Assemble(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                         const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const local_space::ReferenceElement_Data &reference_element_data,
                                         const Polydim::examples::Elastic_PCC_2D::test::I_Test &test) const;

    Performance_Data ComputePerformance(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                        const Gedim::MeshMatricesDAO &mesh,
                                        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                        const local_space::ReferenceElement_Data &reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                         const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const local_space::ReferenceElement_Data &reference_element_data,
                                         const Elastic_PCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Elastic_PCC_2D::test::I_Test &test) const;
};
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
