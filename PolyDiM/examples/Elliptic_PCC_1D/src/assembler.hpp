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
#include "FEM_PCC_1D_LocalSpace.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_1D
{
class Assembler final
{
  public:
    struct Elliptic_PCC_1D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionDirichlet;
    };

    struct PostProcess_Data final
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

  private:
    void ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                           const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                           const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test,
                           Elliptic_PCC_1D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell1DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                         const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test,
                         Elliptic_PCC_1D_Problem_Data &assembler_data) const;

  public:
    Elliptic_PCC_1D_Problem_Data Assemble(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                                          const Gedim::MeshMatricesDAO &mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                          const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                          const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                          const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                          const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace &local_space,
                                          const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                         const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                         const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace &local_space,
                                         const Elliptic_PCC_1D_Problem_Data &assembler_data,
                                         const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test) const;
};
} // namespace Elliptic_PCC_1D
} // namespace examples
} // namespace Polydim

#endif
