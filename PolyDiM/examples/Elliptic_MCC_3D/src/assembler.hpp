#ifndef __assembler_H
#define __assembler_H

#include "Assembler_Utilities.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

#include "DOFsManager.hpp"
#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"

#include "VEM_MCC_3D_LocalSpace_Data.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_3D
{

class Assembler final
{
  public:
    struct Elliptic_MCC_3D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> neumannMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionNeumann;
    };

    struct VEM_Performance_Result final
    {
        struct Cell3D_Performance final
        {
            unsigned int NumBoundaryQuadraturePoints = 0;
            unsigned int NumInternalQuadraturePoints = 0;
            Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data Analysis;
        };

        std::vector<Cell3D_Performance> Cell3DsPerformance;
    };

    struct PostProcess_Data final
    {
        Eigen::VectorXd cell3Ds_numeric_pressure;
        Eigen::VectorXd cell3Ds_exact_pressure;

        Eigen::VectorXd cell3Ds_error_L2_pressure;
        Eigen::VectorXd cell3Ds_super_error_L2_pressure;
        Eigen::VectorXd cell3Ds_norm_L2_pressure;
        double error_L2_pressure;
        double super_error_L2_pressure;
        double norm_L2_pressure;
        Eigen::VectorXd cell3Ds_error_L2_velocity;
        Eigen::VectorXd cell3Ds_norm_L2_velocity;
        double error_L2_velocity;
        double norm_L2_velocity;

        double mesh_size;

        double residual_norm;
    };

  private:
    void ComputeStrongTerm(const unsigned int &cell3DIndex,
                           const Gedim::MeshMatricesDAO &mesh,
                           const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                           const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &local_space_data,
                           const test::I_Test &test,
                           Elliptic_MCC_3D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell3DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                         const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &local_space_data,
                         const test::I_Test &test,
                         Elliptic_MCC_3D_Problem_Data &assembler_data) const;

  public:
    Elliptic_MCC_3D_Problem_Data Assemble(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                                          const Gedim::MeshMatricesDAO &mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                          const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                          const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space,
                                          const Polydim::VEM::MCC::I_VEM_MCC_3D_Pressure_LocalSpace &vem_pressure_space,
                                          const Polydim::examples::Elliptic_MCC_3D::test::I_Test &test) const;

    VEM_Performance_Result ComputeVemPerformance(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                                                 const Gedim::MeshMatricesDAO &mesh,
                                                 const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                                 const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                                 const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                         const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                         const Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                         const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space,
                                         const Polydim::VEM::MCC::I_VEM_MCC_3D_Pressure_LocalSpace &vem_pressure_space,
                                         const Elliptic_MCC_3D_Problem_Data &assembler_data,
                                         const Polydim::examples::Elliptic_MCC_3D::test::I_Test &test) const;
};
} // namespace Elliptic_MCC_3D
} // namespace examples

} // namespace Polydim

#endif
