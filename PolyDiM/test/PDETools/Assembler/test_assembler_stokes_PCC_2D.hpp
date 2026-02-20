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

#ifndef __TEST_assembler_stokes_PCC_2D_H
#define __TEST_assembler_stokes_PCC_2D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Eigen_LUSolver.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "assembler_PCC_2D_functions.hpp"
#include "assembler_PCC_2D_functions_data.hpp"
#include "assembler_PCC_2D_functions_utilities.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
  namespace UnitTesting
  {

    TEST(TEST_assembler_stokes_PCC_2D, TEST_assembler_stokes_PCC_2D_example)
    {
      Gedim::GeometryUtilitiesConfig geometry_utilities_config;
      geometry_utilities_config.Tolerance1D = 1.0e-8;
      geometry_utilities_config.Tolerance2D = 1.0e-12;
      Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

      const std::string exportFolder = "./Export/TEST_assembler_stokes_PCC_2D/TEST_assembler_stokes_PCC_2D_example";
      Gedim::Output::CreateFolder(exportFolder);

      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;
      domain.area = 1.0;
      domain.vertices = Eigen::MatrixXd::Zero(3, 4);
      domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
      domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;
      domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

      Gedim::MeshUtilities mesh_utilities;

      Gedim::MeshMatrices meshData;
      Gedim::MeshMatricesDAO mesh(meshData);

      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometry_utilities,
                                                                  mesh_utilities,
                                                                  Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular,
                                                                  domain,
                                                                  0.001,
                                                                  mesh);

      const double nu = 1.0;
      const unsigned int pressure_method_order = 1;
      const unsigned int velocity_method_order = 2;

      const auto pressure_reference_element_data =
          Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC,
                                                                       pressure_method_order);
      const auto velocity_reference_element_data =
          Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC,
                                                                       velocity_method_order);

      const auto mesh_geometric_data = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(
                                         geometry_utilities,
                                         mesh_utilities,
                                         mesh,
                                         Polydim::PDETools::LocalSpace_PCC_2D::MeshGeometricDataConfigiguration(velocity_reference_element_data));

      Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);
      Polydim::PDETools::DOFs::DOFsManager dofManager;

      std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> pressure_boundary_info = {
        {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

      std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> velocity_boundary_info = {
        {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

      const auto pressure_mesh_dofs_info = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(pressure_reference_element_data, mesh, pressure_boundary_info);
      const auto pressure_dofs_data = dofManager.CreateDOFs_2D(pressure_mesh_dofs_info, mesh_connectivity_data);

      const auto velocity_mesh_dofs_info = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(velocity_reference_element_data, mesh, velocity_boundary_info);
      const auto velocity_dofs_data = dofManager.CreateDOFs_2D(velocity_mesh_dofs_info, mesh_connectivity_data);

      auto pressure_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return sin(2.0 * M_PI * x) * cos(2.0 * M_PI * y);
      };
      auto velocity_x_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return +0.5 * sin(2.0 * M_PI * x) * sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y) * cos(2.0 * M_PI * y);
      };
      auto velocity_y_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return -0.5 * sin(2.0 * M_PI * y) * sin(2.0 * M_PI * y) * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * x);
      };

      auto pressure_gradient_solution_function = [](const double &x, const double &y, const double &z) {
        return std::array<double, 3>({
                                       +2.0 * M_PI * cos(2.0 * M_PI * x) * cos(2.0 * M_PI * y),
                                       -2.0 * M_PI * sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y),
                                       0.0
                                     });
      };

      auto velocity_x_laplacian_solution_function = [](const double &x, const double &y, const double &z) {
        return (+8.0 * M_PI * M_PI * cos(4.0 * M_PI * x) - 4.0 * M_PI * M_PI) * sin(2.0 * M_PI * y) * cos(2.0 * M_PI * y);
      };
      auto velocity_y_laplacian_solution_function = [](const double &x, const double &y, const double &z) {
        return (-8.0 * M_PI * M_PI * cos(4.0 * M_PI * y) + 4.0 * M_PI * M_PI) * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * x);
      };

      auto f_x_function = [&nu, &pressure_gradient_solution_function, &velocity_x_laplacian_solution_function](const double &x, const double &y, const double &z) {
        const auto p_grad = pressure_gradient_solution_function(x, y, z);
        const auto u_x_lap = velocity_x_laplacian_solution_function(x, y, z);

        return - nu * u_x_lap + p_grad.at(0);
      };

      auto f_y_function = [&nu, &pressure_gradient_solution_function, &velocity_y_laplacian_solution_function](const double &x, const double &y, const double &z) {
        const auto p_grad = pressure_gradient_solution_function(x, y, z);
        const auto u_y_lap = velocity_y_laplacian_solution_function(x, y, z);

        return - nu * u_y_lap + p_grad.at(1);
      };

      const auto f_x = PDETools::Assembler_Utilities::PCC_2D::assemble_source_term(geometry_utilities,
                                                                                   mesh,
                                                                                   mesh_geometric_data,
                                                                                   velocity_dofs_data,
                                                                                   velocity_reference_element_data,
                                                                                   velocity_reference_element_data,
                                                                                   f_x_function);

      ASSERT_EQ(velocity_dofs_data.NumberDOFs, f_x.size());

      const auto f_y = PDETools::Assembler_Utilities::PCC_2D::assemble_source_term(geometry_utilities,
                                                                                   mesh,
                                                                                   mesh_geometric_data,
                                                                                   velocity_dofs_data,
                                                                                   velocity_reference_element_data,
                                                                                   velocity_reference_element_data,
                                                                                   f_y_function);
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, f_y.size());

      auto nu_function = [&nu](const double &x, const double &y, const double &z) {
        return nu;
      };

      const auto A = PDETools::Assembler_Utilities::PCC_2D::assemble_diffusion_operator(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        velocity_dofs_data,
                                                                                        velocity_dofs_data,
                                                                                        velocity_reference_element_data,
                                                                                        velocity_reference_element_data,
                                                                                        nu_function);

      ASSERT_EQ(velocity_dofs_data.NumberDOFs, A.operator_dofs.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, A.operator_dofs.size.at(1));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, A.operator_strong.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberStrongs, A.operator_strong.size.at(1));

      auto adv_1_function = [](const double &x, const double &y, const double &z) {
        return std::array<double, 3>({ 1.0, 0.0, 0.0 });
      };
      auto adv_2_function = [](const double &x, const double &y, const double &z) {
        return std::array<double, 3>({ 0.0, 1.0, 0.0 });
      };

      const auto B_x = PDETools::Assembler_Utilities::PCC_2D::assemble_advection_operator(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          pressure_dofs_data,
                                                                                          velocity_dofs_data,
                                                                                          pressure_reference_element_data,
                                                                                          velocity_reference_element_data,
                                                                                          adv_1_function);
      const auto B_y = PDETools::Assembler_Utilities::PCC_2D::assemble_advection_operator(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          pressure_dofs_data,
                                                                                          velocity_dofs_data,
                                                                                          pressure_reference_element_data,
                                                                                          velocity_reference_element_data,
                                                                                          adv_2_function);

      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_x.operator_dofs.size.at(0));
      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_x.operator_dofs.size.at(1));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_x.operator_strong.size.at(0));
      ASSERT_EQ(pressure_dofs_data.NumberStrongs, B_x.operator_strong.size.at(1));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_y.operator_dofs.size.at(0));
      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_y.operator_dofs.size.at(1));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_y.operator_strong.size.at(0));
      ASSERT_EQ(pressure_dofs_data.NumberStrongs, B_y.operator_strong.size.at(1));

      auto pressure_strong_function =
          [&pressure_exact_solution_function](const unsigned int marker, const double &x, const double &y, const double &z) {
        if (marker != 1)
          throw std::runtime_error("marker not managed");

        return pressure_exact_solution_function(x, y, z);
      };

      const auto p_strong = PDETools::Assembler_Utilities::PCC_2D::assemble_strong_solution(geometry_utilities,
                                                                                            mesh,
                                                                                            mesh_geometric_data,
                                                                                            pressure_mesh_dofs_info,
                                                                                            pressure_dofs_data,
                                                                                            pressure_reference_element_data,
                                                                                            pressure_strong_function);

      ASSERT_EQ(pressure_dofs_data.NumberStrongs, p_strong.size());

      Eigen::VectorXd u_x_strong = Eigen::VectorXd::Zero(velocity_dofs_data.NumberStrongs);
      Eigen::VectorXd u_y_strong = Eigen::VectorXd::Zero(velocity_dofs_data.NumberStrongs);

      {
        const unsigned int tot_dofs = 2 * velocity_dofs_data.NumberDOFs + pressure_dofs_data.NumberDOFs;
        const unsigned int tot_strongs = 2 * velocity_dofs_data.NumberStrongs + pressure_dofs_data.NumberStrongs;

        Eigen::VectorXd p_numeric;
        Eigen::VectorXd u_x_numeric;
        Eigen::VectorXd u_y_numeric;

        {
          Eigen::VectorXd f_concat = Eigen::VectorXd::Zero(tot_dofs);
          f_concat.segment(0, velocity_dofs_data.NumberDOFs) = f_x;
          f_concat.segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs) = f_y;

          const auto f = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(f_concat);

          Eigen::VectorXd u_D_concat = Eigen::VectorXd::Zero(tot_strongs);
          u_D_concat.segment(0, velocity_dofs_data.NumberStrongs) = u_x_strong;
          u_D_concat.segment(velocity_dofs_data.NumberStrongs, velocity_dofs_data.NumberStrongs) = u_y_strong;
          u_D_concat.segment(2 * velocity_dofs_data.NumberStrongs, pressure_dofs_data.NumberStrongs) = p_strong;
          const auto u_D = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(u_D_concat);

          const auto J_A_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                               PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(A.operator_dofs,
                                                                                           { tot_dofs, tot_dofs },
                                                                                           { 0, 0 }));
          const auto J_A_x_D = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                 PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(A.operator_strong,
                                                                                             { tot_dofs, tot_strongs },
                                                                                             { 0, 0 }));
          const auto J_A_y = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                               PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(A.operator_dofs,
                                                                                           { tot_dofs, tot_dofs },
                                                                                           { velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs }));
          const auto J_A_y_D = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                 PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(A.operator_strong,
                                                                                             { tot_dofs, tot_strongs },
                                                                                             { velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberStrongs }));
          const auto J_B_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                               PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_x.operator_dofs,
                                                                                           { tot_dofs, tot_dofs },
                                                                                           { 0, 2 * velocity_dofs_data.NumberDOFs }));
          const auto J_B_x_D = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                 PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_x.operator_strong,
                                                                                             { tot_dofs, tot_strongs },
                                                                                             { 0, 2 * velocity_dofs_data.NumberStrongs }));
          const auto J_B_y = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                               PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_y.operator_dofs,
                                                                                           { tot_dofs, tot_dofs },
                                                                                           { velocity_dofs_data.NumberDOFs, 2 * velocity_dofs_data.NumberDOFs }));
          const auto J_B_y_D = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                 PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_y.operator_strong,
                                                                                             { tot_dofs, tot_strongs },
                                                                                             { velocity_dofs_data.NumberDOFs, 2 * velocity_dofs_data.NumberStrongs }));
          const auto J_BT_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_x.operator_dofs,
                                                                                            { tot_dofs, tot_dofs },
                                                                                            { 0, 2 * velocity_dofs_data.NumberDOFs },
                                                                                            true));
          const auto J_BT_y = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                                PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_y.operator_dofs,
                                                                                            { tot_dofs, tot_dofs },
                                                                                            { velocity_dofs_data.NumberDOFs, 2 * velocity_dofs_data.NumberDOFs },
                                                                                            true));

          const auto J = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(J_A_x + J_A_y - J_B_x - J_B_y - J_BT_x - J_BT_y);
          const auto J_D = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(J_A_x_D + J_A_y_D - J_B_x_D - J_B_y_D);

          auto rhs = f;
          rhs.SubtractionMultiplication(J_D, u_D);

          Gedim::Eigen_Array<> u;
          u.SetSize(tot_dofs);

          Gedim::Eigen_LUSolver solver;
          solver.Initialize(J);
          solver.Solve(rhs, u);

          const auto numeric_solution = PDETools::Assembler_Utilities::PCC_2D::to_VectorXd(u);

          u_x_numeric = u.Segment(0, velocity_dofs_data.NumberDOFs);
          u_y_numeric = u.Segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs);
          p_numeric = u.Segment(2 * velocity_dofs_data.NumberDOFs, pressure_dofs_data.NumberDOFs);

        }

        const auto u_x_on_cell0Ds = PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh,
                                                                                                       velocity_dofs_data,
                                                                                                       u_x_numeric,
                                                                                                       u_x_strong,
                                                                                                       velocity_x_exact_solution_function);
        const auto u_y_on_cell0Ds = PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh,
                                                                                                       velocity_dofs_data,
                                                                                                       u_y_numeric,
                                                                                                       u_y_strong,
                                                                                                       velocity_y_exact_solution_function);
        const auto p_on_cell0Ds = PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh,
                                                                                                     pressure_dofs_data,
                                                                                                     p_numeric,
                                                                                                     p_strong,
                                                                                                     pressure_exact_solution_function);


        const auto p_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                      mesh,
                                                                                      mesh_geometric_data,
                                                                                      pressure_dofs_data,
                                                                                      pressure_reference_element_data,
                                                                                      p_numeric,
                                                                                      p_strong,
                                                                                      pressure_exact_solution_function);

        const auto u_x_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                      mesh,
                                                                                      mesh_geometric_data,
                                                                                      velocity_dofs_data,
                                                                                      velocity_reference_element_data,
                                                                                      u_x_numeric,
                                                                                      u_x_strong,
                                                                                      velocity_x_exact_solution_function);
        const auto u_y_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                      mesh,
                                                                                      mesh_geometric_data,
                                                                                      velocity_dofs_data,
                                                                                      velocity_reference_element_data,
                                                                                      u_y_numeric,
                                                                                      u_y_strong,
                                                                                      velocity_y_exact_solution_function);

        //        const auto u_on_quadrature = PDETools::Assembler_Utilities::PCC_2D::evaluate_solution_on_quadrature_points(geometry_utilities,
        //                                                                                                                   mesh,
        //                                                                                                                   mesh_geometric_data,
        //                                                                                                                   trial_dofs_data,
        //                                                                                                                   trial_reference_element_data,
        //                                                                                                                   numeric_solution,
        //                                                                                                                   strong_solution,
        //                                                                                                                   exact_solution_function,
        //                                                                                                                   exact_gradient_solution_function);

                {
                  Gedim::VTKUtilities exporter;

                  exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                                       mesh.Cell2DsVertices(),
                                       {
                                         {
                                           "u_x_numeric",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(u_x_on_cell0Ds.numeric_solution.size()),
                                           u_x_on_cell0Ds.numeric_solution.data()
                                         },
                                         {
                                           "u_x_exact",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(u_x_on_cell0Ds.exact_solution.size()),
                                           u_x_on_cell0Ds.exact_solution.data()
                                         },
                                         {
                                           "u_y_numeric",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(u_y_on_cell0Ds.numeric_solution.size()),
                                           u_y_on_cell0Ds.numeric_solution.data()
                                         },
                                         {
                                           "u_y_exact",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(u_y_on_cell0Ds.exact_solution.size()),
                                           u_y_on_cell0Ds.exact_solution.data()
                                         },
                                         {
                                           "p_numeric",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(p_on_cell0Ds.numeric_solution.size()),
                                           p_on_cell0Ds.numeric_solution.data()
                                         },
                                         {
                                           "p_exact",
                                           Gedim::VTPProperty::Formats::Points,
                                           static_cast<unsigned int>(p_on_cell0Ds.exact_solution.size()),
                                           p_on_cell0Ds.exact_solution.data()
                                         },
                                         {
                                           "u_x_error_L2",
                                           Gedim::VTPProperty::Formats::Cells,
                                           static_cast<unsigned int>(u_x_error_L2.cell2Ds_error_L2.size()),
                                           u_x_error_L2.cell2Ds_error_L2.data()
                                         },
                                         {
                                           "u_y_error_L2",
                                           Gedim::VTPProperty::Formats::Cells,
                                           static_cast<unsigned int>(u_y_error_L2.cell2Ds_error_L2.size()),
                                           u_y_error_L2.cell2Ds_error_L2.data()
                                         },
                                         {
                                           "p_error_L2",
                                           Gedim::VTPProperty::Formats::Cells,
                                           static_cast<unsigned int>(p_error_L2.cell2Ds_error_L2.size()),
                                           p_error_L2.cell2Ds_error_L2.data()
                                         }
                                       });
                  exporter.Export(exportFolder + "/solution.vtu");
                }

        //        {
        //          Gedim::VTKUtilities exporter;

        //          exporter.AddPoints(u_on_quadrature.quadrature_points,
        //                             {
        //                               {
        //                                 "numeric_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.numeric_solution.size()),
        //                                 u_on_quadrature.numeric_solution.data()
        //                               },
        //                               {
        //                                 "numeric_x_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.numeric_gradient_solution.at(0).size()),
        //                                 u_on_quadrature.numeric_gradient_solution.at(0).data()
        //                               },
        //                               {
        //                                 "numeric_y_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.numeric_gradient_solution.at(1).size()),
        //                                 u_on_quadrature.numeric_gradient_solution.at(1).data()
        //                               },
        //                               {
        //                                 "exact_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.exact_solution.size()),
        //                                 u_on_quadrature.exact_solution.data()
        //                               },
        //                               {
        //                                 "exact_x_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.exact_gradient_solution.at(0).size()),
        //                                 u_on_quadrature.exact_gradient_solution.at(0).data()
        //                               },
        //                               {
        //                                 "exact_y_solution",
        //                                 Gedim::VTPProperty::Formats::Points,
        //                                 static_cast<unsigned int>(u_on_quadrature.exact_gradient_solution.at(1).size()),
        //                                 u_on_quadrature.exact_gradient_solution.at(1).data()
        //                               },
        //                             });

        //          exporter.Export(exportFolder + "/solution_on_quadrature.vtu");
        //        }

        //        //std::cout.precision(2);
        //        //std::cout<< std::scientific<< "A: "<< A<< std::endl;
        //        //std::cout<< std::scientific<< "A_D: "<< A_D<< std::endl;
        //        //std::cout << std::scientific << "f: " << f << std::endl;
        //        //std::cout << std::scientific << "w_t: " << w_t << std::endl;
        //        //std::cout << std::scientific << "r: " << rhs << std::endl;
        //        //std::cout << std::scientific << "u: " << u << std::endl;
        //        //std::cout << std::scientific << "u_ex: " << exact_solution.exact_solution.transpose() << std::endl;
        //        //std::cout << std::scientific << "u_D: " << u_D << std::endl;
        //        //std::cout << std::scientific << "u_ex_D: " << exact_solution.exact_solution_strong.transpose() << std::endl;
        //        //std::cout << std::scientific << "err_L2: " << (error_L2.error_L2 / error_L2.exact_norm_L2) << std::endl;
        //        //std::cout << std::scientific << "err_H1: " << (error_H1.error_H1 / error_H1.exact_norm_H1) << std::endl;

        std::cout.precision(2);
        std::cout<< std::scientific<< "u_x_errorL2: "<< (u_x_error_L2.error_L2 / u_x_error_L2.numeric_norm_L2)<< " ";
        std::cout<< std::scientific<< "u_y_errorL2: "<< (u_y_error_L2.error_L2 / u_y_error_L2.numeric_norm_L2)<< " ";
        std::cout<< std::scientific<< "p_errorL2: "<< (p_error_L2.error_L2 / p_error_L2.numeric_norm_L2)<< std::endl;

        //        ASSERT_TRUE(error_L2.error_L2 < 1.0e-13 * error_L2.exact_norm_L2);
        //        ASSERT_TRUE(error_H1.error_H1 < 1.0e-13 * error_H1.exact_norm_H1);
      }
    }

  } // namespace UnitTesting
} // namespace Polydim

#endif
