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

#include "DOFsManager.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
class Assembler final
{
  public:
    struct Parabolic_PCC_2D_Initial_Data final
    {
        Gedim::Eigen_Array<> initial_condition;
        Gedim::Eigen_Array<> initial_condition_dirichlet;
    };

    struct Parabolic_PCC_2D_Static_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_SparseArray<> globalMatrixM;
        Gedim::Eigen_SparseArray<> dirichletMatrixM;
    };

    struct Parabolic_PCC_2D_Problem_Data final
    {
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionDirichlet;
    };

    struct Performance_Data final
    {
        std::vector<Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data> Cell2DsPerformance;
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
    void ComputeStrongTerm(const unsigned int cell2D_index,
                           const Gedim::MeshMatricesDAO &mesh,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                           const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                           const test::I_Test &test,
                           const double &time_value,
                           Parabolic_PCC_2D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                         const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                         const double &time_value,
                         Parabolic_PCC_2D_Problem_Data &assembler_data) const;

  public:
    Parabolic_PCC_2D_Static_Problem_Data StaticAssemble(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                                         const Gedim::MeshMatricesDAO &mesh,
                                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test) const;

    Parabolic_PCC_2D_Problem_Data Assemble(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                           const Gedim::MeshMatricesDAO &mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                           const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                                           const Assembler::Parabolic_PCC_2D_Static_Problem_Data& static_assembler_data,
                                           const double &time_value) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                         const Parabolic_PCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                                         const Gedim::Eigen_SparseArray<>& A,
                                         const Gedim::Eigen_Array<>& rhs,
                                         const double &time_value) const;

    Parabolic_PCC_2D_Initial_Data ComputeInitalCondition(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                const Gedim::IMeshDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test) const;
};
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
