#ifndef __assembler_H
#define __assembler_H

#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"
#include "DOFsManager.hpp"

#include "VEM_MCC_2D_Velocity_LocalSpace_Data.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{

class Assembler final
{
public:

    struct Elliptic_MCC_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> neumannMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionNeumann;
    };

    struct VEM_Performance_Result final
    {
        struct Cell2D_Performance final
        {
            unsigned int NumBoundaryQuadraturePoints = 0;
            unsigned int NumInternalQuadraturePoints = 0;
            Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data Analysis;
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
        Eigen::VectorXd cell2Ds_error_L2_velocity;
        Eigen::VectorXd cell2Ds_norm_L2_velocity;
        double error_L2_velocity;
        double norm_L2_velocity;

        double mesh_size;

        double residual_norm;
    };

private:

    void ComputeStrongTerm(const Gedim::MeshMatricesDAO& mesh,
                           const unsigned int &cell2DIndex,
                           const std::vector<bool> &cell2DEdgeDirections,
                           const Eigen::MatrixXd &boundaryQuadraturePoints,
                           const Eigen::VectorXd &boundaryQuadratureWeights,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                           const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                           const test::I_Test &test,
                           Elliptic_MCC_2D_Problem_Data& assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO& mesh,
                         const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry& polygon,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                         const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &local_space_data,
                         const test::I_Test &test,
                         Elliptic_MCC_2D_Problem_Data& assembler_data) const;

public:
    Elliptic_MCC_2D_Problem_Data Assemble(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                                          const Gedim::MeshMatricesDAO& mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo>& mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData>& dofs_data,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                          const Polydim::examples::Elliptic_MCC_2D::test::I_Test& test) const;

    VEM_Performance_Result ComputeVemPerformance(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                                                 const Gedim::MeshMatricesDAO& mesh,
                                                 const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                 const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                                         const Gedim::MeshMatricesDAO& mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData>& dofs_data,
                                         const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                         const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                         const Elliptic_MCC_2D_Problem_Data& assembler_data,
                                         const Polydim::examples::Elliptic_MCC_2D::test::I_Test& test) const;
};
}
}

}

#endif
