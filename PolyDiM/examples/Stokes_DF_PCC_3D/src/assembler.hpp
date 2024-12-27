#ifndef __assembler_H
#define __assembler_H

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

#include "DOFsManager.hpp"
#include "VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_PerformanceAnalysis.hpp"

#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
class Assembler final
{
  public:
    struct Stokes_DF_PCC_3D_Problem_Data final
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
            Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis_Data Analysis;
        };

        std::vector<Cell2D_Performance> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        std::vector<Eigen::VectorXd> cell0Ds_numeric_pressure;
        std::vector<Eigen::VectorXd> cell0Ds_exact_pressure;

        Eigen::VectorXd cell2Ds_error_L2_pressure;
        Eigen::VectorXd cell2Ds_norm_L2_pressure;
        double error_L2_pressure;
        double norm_L2_pressure;
        Eigen::VectorXd cell2Ds_error_H1_velocity;
        Eigen::VectorXd cell2Ds_norm_H1_velocity;
        double error_H1_velocity;
        double norm_H1_velocity;

        double mesh_size;

        double residual_norm;
    };

  private:
    void ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                           const std::vector<size_t> &offsetStrongs,
                           const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                           const test::I_Test &test,
                           Stokes_DF_PCC_3D_Problem_Data &assembler_data) const;

  public:
    Stokes_DF_PCC_3D_Problem_Data Assemble(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                                           const Gedim::MeshMatricesDAO &mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                           const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                           const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                           const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                           const Polydim::examples::Stokes_DF_PCC_3D::test::I_Test &test) const;

    VEM_Performance_Result ComputeVemPerformance(
        const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                         const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                         const Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                                         const Polydim::examples::Stokes_DF_PCC_3D::test::I_Test &test) const;
};
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
