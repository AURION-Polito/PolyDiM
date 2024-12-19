#ifndef __assembler_H
#define __assembler_H

#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "DOFsManager.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"


namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
class Assembler final
{
public:
    struct Elliptic_PCC_3D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionDirichlet;
    };

    struct VEM_Performance_Result final
    {
        struct Cell3D_Performance final
        {
            unsigned int NumBoundaryQuadraturePoints = 0;
            unsigned int NumInternalQuadraturePoints = 0;
            Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
        };

        std::vector<Cell3D_Performance> Cell3DsPerformance;
    };

    struct PostProcess_Data final
    {
        Eigen::VectorXd cell0Ds_numeric;
        Eigen::VectorXd cell0Ds_exact;

        Eigen::VectorXd cell3Ds_error_L2;
        Eigen::VectorXd cell3Ds_norm_L2;
        double error_L2;
        double norm_L2;
        Eigen::VectorXd cell3Ds_error_H1;
        Eigen::VectorXd cell3Ds_norm_H1;
        double error_H1;
        double norm_H1;

        double mesh_size;

        double residual_norm;
    };

private:
    void ComputeStrongTerm(const unsigned int &cell3DIndex,
                           const Gedim::MeshMatricesDAO& mesh,
                           const std::vector<VEM::PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                           const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                           const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                           const Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Data &local_space_data,
                           const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test,
                           Elliptic_PCC_3D_Problem_Data& assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell3DIndex,
                         const Gedim::MeshMatricesDAO& mesh,
                         const Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                         const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                         const VEM::PCC::VEM_PCC_3D_LocalSpace_Data &local_space_data,
                         const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test,
                         Elliptic_PCC_3D_Problem_Data& assembler_data) const;

public:
    Elliptic_PCC_3D_Problem_Data Assemble(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                          const Gedim::MeshMatricesDAO& mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                          const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                          const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                          const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                          const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                                          const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test) const;

    VEM_Performance_Result ComputeVemPerformance(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                                 const Gedim::MeshMatricesDAO& mesh,
                                                 const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                                 const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                                 const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                         const Gedim::MeshMatricesDAO& mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                         const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                         const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                                         const Elliptic_PCC_3D_Problem_Data& assembler_data,
                                         const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test) const;

};
}
}
}

#endif
