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
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "local_space.hpp"
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

    struct Performance_Data final
    {
        std::vector<local_space::Performance_Data> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        Eigen::VectorXd cell2Ds_numeric_pressure;
        Eigen::VectorXd cell2Ds_exact_pressure;

        Eigen::VectorXd cell2Ds_error_L2_pressure;
        Eigen::VectorXd cell2Ds_super_error_L2_pressure;
        Eigen::VectorXd cell2Ds_norm_L2_pressure;
        double error_L2_pressure;
        double super_error_L2_pressure;
        double norm_L2_pressure;
        Eigen::VectorXd cell2Ds_error_L2_velocity;
        Eigen::VectorXd cell2Ds_norm_L2_velocity;
        double error_L2_velocity;
        double norm_L2_velocity;

        double mesh_size;

        double residual_norm;
    };

  private:
    void ComputeStrongTerm(const unsigned int &cell2DIndex,
                           const Gedim::MeshMatricesDAO &mesh,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                           const local_space::ReferenceElement_Data &reference_element_data,
                           const local_space::LocalSpace_Data &local_space_data,
                           const test::I_Test &test,
                           Elliptic_MCC_2D_Problem_Data &assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                         const local_space::ReferenceElement_Data &reference_element_data,
                         const local_space::LocalSpace_Data &local_space_data,
                         const test::I_Test &test,
                         Elliptic_MCC_2D_Problem_Data &assembler_data) const;

  public:
    Elliptic_MCC_2D_Problem_Data Assemble(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                          const Gedim::MeshMatricesDAO &mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                          const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                          const local_space::ReferenceElement_Data &reference_element_data,
                                          const Polydim::examples::Elliptic_MCC_2D::test::I_Test &test) const;

    Assembler::Performance_Data ComputePerformance(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                                   const Gedim::MeshMatricesDAO &mesh,
                                                   const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                   const local_space::ReferenceElement_Data &reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                         const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                         const local_space::ReferenceElement_Data &reference_element_data,
                                         const Elliptic_MCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Elliptic_MCC_2D::test::I_Test &test) const;
};
} // namespace Elliptic_MCC_2D
} // namespace examples

} // namespace Polydim

#endif
