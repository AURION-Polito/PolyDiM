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

#ifndef __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_HPP
#define __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_HPP

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshUtilities.hpp"

#include "assembler_PCC_2D_functions_data.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
namespace PCC_2D
{
// ***************************************************************************
Eigen::VectorXd assembler_source_term(
    const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<double(const double &, const double &, const double &, const Eigen::VectorXd &)> source_term_function);
// ***************************************************************************
Variational_Operator assembler_elliptic_operator(const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<double(const double &, const double &, const double &)> diffusion_term_function,
    const std::function<std::array<double, 3> (const double&, const double&, const double&)> advection_term_function,
    const std::function<double (const double&, const double&, const double&)> reaction_term_function);
// ***************************************************************************
inline Variational_Operator assembler_diffusion_operator(const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<double(const double &, const double &, const double &)> diffusion_term_function)
{
  return assembler_elliptic_operator(geometry_utilities,
                                     mesh,
                                     mesh_geometric_data,
                                     dofs_data,
                                     reference_element_data,
                                     diffusion_term_function,
                                     nullptr,
                                     nullptr);
}
// ***************************************************************************
inline Variational_Operator assembler_reaction_operator(const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<double (const double&, const double&, const double&)> reaction_term_function)
{
  return assembler_elliptic_operator(geometry_utilities,
                                     mesh,
                                     mesh_geometric_data,
                                     dofs_data,
                                     reference_element_data,
                                     nullptr,
                                     nullptr,
                                     reaction_term_function);
}
// ***************************************************************************
inline Variational_Operator assembler_advection_operator(const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<std::array<double, 3> (const double&, const double&, const double&)> advection_term_function)
{
  return assembler_elliptic_operator(geometry_utilities,
                                     mesh,
                                     mesh_geometric_data,
                                     dofs_data,
                                     reference_element_data,
                                     nullptr,
                                     advection_term_function,
                                     nullptr);
}
// ***************************************************************************
Eigen::VectorXd assembler_strong_solution(
    const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const std::function<double(const unsigned int, const double &, const double &, const double &)> strong_solution_function);
// ***************************************************************************
Exact_Solution_Data assembler_exact_solution(const Gedim::GeometryUtilities &geometry_utilities,
                                             const Gedim::MeshMatricesDAO &mesh,
                                             const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                             const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                             const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                             const std::function<double(const double &, const double &, const double &)> exact_solution_function);
// ***************************************************************************
Eigen::VectorXd assembler_weak_term(const Gedim::GeometryUtilities &geometry_utilities,
                                    const Gedim::MeshMatricesDAO &mesh,
                                    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                    const DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                    const std::function<double(const unsigned int, const double &, const double &, const double &)> weak_term_function);
// ***************************************************************************
Post_Process_Data_Cell0Ds assembler_extract_cell0Ds(const Gedim::MeshMatricesDAO &mesh,
                                                    const DOFs::DOFsManager::DOFsData &dofs_data,
                                                    const Eigen::VectorXd &numerical_solution,
                                                    const Eigen::VectorXd &numerical_solution_strong,
                                                    const std::function<double(const double &, const double &, const double &)> exact_solution_function);
// ***************************************************************************
Post_Process_Data_ErrorL2 assembler_error_L2(const Gedim::GeometryUtilities &geometry_utilities,
                                             const Gedim::MeshMatricesDAO &mesh,
                                             const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                             const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                             const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Eigen::VectorXd &numerical_solution,
                                             const Eigen::VectorXd &numerical_solution_strong,
                                             const std::function<double(const double &, const double &, const double &)> exact_solution_function);
// ***************************************************************************
Post_Process_Data_ErrorH1 assembler_error_H1(const Gedim::GeometryUtilities &geometry_utilities,
                                             const Gedim::MeshMatricesDAO &mesh,
                                             const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                             const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                             const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Eigen::VectorXd &numerical_solution,
                                             const Eigen::VectorXd &numerical_solution_strong,
                                             const std::function<std::array<double, 3>(const double &, const double &, const double &)> exact_gradient_solution_function);
// ***************************************************************************
} // namespace PCC_2D
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
