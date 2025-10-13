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

#ifndef __assembler_H
#define __assembler_H

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "EllipticEquation.hpp"
#include "FEM_PCC_1D_LocalSpace.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_BulkFace_2D
{
class Assembler final
{

    Polydim::PDETools::Equations::EllipticEquation equation;

  public:
    struct Elliptic_PCC_BF_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> initial_solution;
    };

    struct Performance_Data_2D final
    {
        std::vector<Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data> Cell2DsPerformance;
    };

    struct PostProcess_Data_2D final
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
    };

    struct PostProcess_Data_1D final
    {
        Eigen::VectorXd cell0Ds_numeric;
        Eigen::VectorXd cell0Ds_exact;

        Eigen::VectorXd cell1Ds_error_L2;
        Eigen::VectorXd cell1Ds_norm_L2;
        double error_L2;
        double norm_L2;
        Eigen::VectorXd cell1Ds_error_H1;
        Eigen::VectorXd cell1Ds_norm_H1;
        double error_H1;
        double norm_H1;

        double mesh_size;

        double residual_norm;
    };

    struct PostProcess_Data final
    {
        PostProcess_Data_2D post_process_data_2D;
        PostProcess_Data_1D post_process_data_1D;
        double residual_norm;
    };

  private:
    void Assemble_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                     const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                     Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const;

    void Assemble_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                     const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                     Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const;

    void PostProcessSolution_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                                Assembler::PostProcess_Data_2D &result) const;

    void PostProcessSolution_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                                Assembler::PostProcess_Data_1D &result) const;

    void ComputeTransitionMatrices(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                   const Gedim::MeshMatricesDAO &mesh_1D,
                                   const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                                   const Gedim::MeshUtilities::ExtractMeshData &extract_data,
                                   const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                   const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                   const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                   const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                                   const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                                   Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const;

    void ComputeInitalCondition_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                   const Gedim::IMeshDAO &mesh,
                                   const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                   const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                   const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                   const test::I_Test &test,
                                   Gedim::Eigen_Array<> &initial_solution) const;

    void ComputeInitalCondition_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                   const Gedim::IMeshDAO &mesh,
                                   const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                   const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                   const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                   const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                   const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                   const test::I_Test &test,
                                   Gedim::Eigen_Array<> &initial_condition) const;

  public:
    Elliptic_PCC_BF_2D_Problem_Data Solve(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                          const std::vector<double> &time_steps,
                                          const Gedim::MeshMatricesDAO &mesh_2D,
                                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
                                          const Gedim::MeshMatricesDAO &mesh_1D,
                                          const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                                          const Gedim::MeshUtilities::ExtractMeshData &extract_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                          const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                          const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
                                          const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                                          const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test) const;

    Performance_Data_2D ComputePerformance_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                              const Gedim::MeshMatricesDAO &mesh,
                                              const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                              const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data) const;

    Assembler::PostProcess_Data PostProcessSolution(
        const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh_2D,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
        const Gedim::MeshMatricesDAO &mesh_1D,
        const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
        const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
        const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
        const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
        const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
        const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
        const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
        const Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const;
};
} // namespace Elliptic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim

#endif
