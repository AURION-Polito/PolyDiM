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

#ifndef __TEST_assembler_navier_stokes_PCC_2D_H
#define __TEST_assembler_navier_stokes_PCC_2D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Eigen_LUSolver.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler_PCC_2D_functions.hpp"
#include "assembler_PCC_2D_functions_data.hpp"
#include "assembler_PCC_2D_functions_utilities.hpp"

namespace Polydim
{
  namespace UnitTesting
  {

    TEST(TEST_assembler_navier_stokes_PCC_2D, TEST_assembler_navier_stokes_PCC_2D_example)
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
                                                                  0.01,
                                                                  mesh);

      const double nu = 0.1;
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

      const auto pressure_mesh_dofs_info =
          Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(pressure_reference_element_data, mesh, pressure_boundary_info);
      const auto pressure_dofs_data = dofManager.CreateDOFs_2D(pressure_mesh_dofs_info, mesh_connectivity_data);

      const auto velocity_mesh_dofs_info =
          Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(velocity_reference_element_data, mesh, velocity_boundary_info);
      const auto velocity_dofs_data = dofManager.CreateDOFs_2D(velocity_mesh_dofs_info, mesh_connectivity_data);

      auto pressure_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return std::numbers::pi * std::numbers::pi * sin(2.0 * std::numbers::pi * x) *
            cos(2.0 * std::numbers::pi * y);
      };
      auto velocity_x_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return 0.5 * sin(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * x) *
            sin(2.0 * std::numbers::pi * y) * cos(2.0 * std::numbers::pi * y);
      };
      auto velocity_y_exact_solution_function = [](const double &x, const double &y, const double &z) {
        return -0.5 * sin(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * y) *
            sin(2.0 * std::numbers::pi * x) * cos(2.0 * std::numbers::pi * x);
      };

      auto pressure_gradient_solution_function = [](const double &x, const double &y, const double &z) {
        return std::array<double, 3>(
              {
                2.0 * std::numbers::pi * std::numbers::pi * std::numbers::pi *
                cos(2.0 * std::numbers::pi * x) *
                cos(2.0 * std::numbers::pi * y),
                -2.0 * std::numbers::pi * std::numbers::pi * std::numbers::pi *
                sin(2.0 * std::numbers::pi * x) *
                sin(2.0 * std::numbers::pi * y),
                0.0
              });
      };

      auto velocity_x_laplacian_solution_function = [](const double &x, const double &y, const double &z) {
        return 2.0 * std::numbers::pi * std::numbers::pi * sin(4.0 * std::numbers::pi * y) *
            (2.0 * cos(4.0 * std::numbers::pi * x) - 1.0);
      };
      auto velocity_y_laplacian_solution_function = [](const double &x, const double &y, const double &z) {
        return -2.0 * std::numbers::pi * std::numbers::pi * sin(4.0 * std::numbers::pi * x) *
            (2.0 * cos(4.0 * std::numbers::pi * y) - 1.0);
      };

      auto f_x_function = [&nu,
                          &pressure_gradient_solution_function,
                          &velocity_x_laplacian_solution_function](const double &x, const double &y, const double &z) {
        const auto p_grad = pressure_gradient_solution_function(x, y, z);
        const auto u_x_lap = velocity_x_laplacian_solution_function(x, y, z);

        const double c_x = (std::numbers::pi * cos(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * x) *
                            sin(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * x) *
                            sin(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * y)) *
                           0.5;

        return -nu * u_x_lap + c_x + p_grad.at(0);
      };

      auto f_y_function = [&nu,
                          &pressure_gradient_solution_function,
                          &velocity_y_laplacian_solution_function](const double &x, const double &y, const double &z) {
        const auto p_grad = pressure_gradient_solution_function(x, y, z);
        const auto u_y_lap = velocity_y_laplacian_solution_function(x, y, z);

        const double c_y = (std::numbers::pi * cos(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * x) *
                            sin(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * y) *
                            sin(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * y)) *
                           0.5;

        return -nu * u_y_lap + c_y + p_grad.at(1);
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

      auto nu_function = [&nu](const double &x, const double &y, const double &z) { return nu; };

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
        return std::array<double, 3>({1.0, 0.0, 0.0});
      };
      auto adv_2_function = [](const double &x, const double &y, const double &z) {
        return std::array<double, 3>({0.0, 1.0, 0.0});
      };

      const auto B_x = PDETools::Assembler_Utilities::PCC_2D::assemble_advection_operator(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          velocity_dofs_data,
                                                                                          pressure_dofs_data,
                                                                                          velocity_reference_element_data,
                                                                                          pressure_reference_element_data,
                                                                                          adv_1_function);
      const auto B_y = PDETools::Assembler_Utilities::PCC_2D::assemble_advection_operator(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          velocity_dofs_data,
                                                                                          pressure_dofs_data,
                                                                                          velocity_reference_element_data,
                                                                                          pressure_reference_element_data,
                                                                                          adv_2_function);

      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_x.operator_dofs.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_x.operator_dofs.size.at(1));
      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_x.operator_strong.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberStrongs, B_x.operator_strong.size.at(1));
      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_y.operator_dofs.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberDOFs, B_y.operator_dofs.size.at(1));
      ASSERT_EQ(pressure_dofs_data.NumberDOFs, B_y.operator_strong.size.at(0));
      ASSERT_EQ(velocity_dofs_data.NumberStrongs, B_y.operator_strong.size.at(1));

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

      double residual_norm = 1.0;
      double solution_norm = 1.0;
      double newton_tol = 1.0e-6;
      double max_iterations = 7;
      double num_iteration = 1;

      const unsigned int tot_dofs = 2 * velocity_dofs_data.NumberDOFs + pressure_dofs_data.NumberDOFs;
      const unsigned int tot_strongs = 2 * velocity_dofs_data.NumberStrongs + pressure_dofs_data.NumberStrongs;

      Eigen::VectorXd u_D_concat = Eigen::VectorXd::Zero(tot_strongs);
      u_D_concat.segment(0, velocity_dofs_data.NumberStrongs) = u_x_strong;
      u_D_concat.segment(velocity_dofs_data.NumberStrongs, velocity_dofs_data.NumberStrongs) = u_y_strong;
      u_D_concat.segment(2 * velocity_dofs_data.NumberStrongs, pressure_dofs_data.NumberStrongs) = p_strong;

      Eigen::SparseMatrix<double> J_S;
      Eigen::VectorXd f_S;
      Gedim::Eigen_SparseArray<> J_Stokes;
      {
        const auto J_A_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                             PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(A.operator_dofs, {tot_dofs, tot_dofs}, {0, 0}));
        const auto J_A_y = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(
                                                                                    A.operator_dofs,
                                                                                    {tot_dofs, tot_dofs},
                                                                                    {velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs}));
        const auto J_B_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                             PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_x.operator_dofs,
                                                                                         {tot_dofs, tot_dofs},
                                                                                         {2 * velocity_dofs_data.NumberDOFs, 0}));
        const auto J_B_y = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(
                                                                                    B_y.operator_dofs,
                                                                                    {tot_dofs, tot_dofs},
                                                                                    {2 * velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs}));
        const auto J_BT_x = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                              PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(B_x.operator_dofs,
                                                                                          {tot_dofs, tot_dofs},
                                                                                          {2 * velocity_dofs_data.NumberDOFs, 0},
                                                                                          true));
        const auto J_BT_y =
            PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(
                                                                     B_y.operator_dofs,
                                                                     {tot_dofs, tot_dofs},
                                                                     {2 * velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs},
                                                                     true));

        J_S = J_A_x + J_A_y - J_B_x - J_B_y - J_BT_x - J_BT_y;
        J_Stokes = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(J_S);

        f_S = Eigen::VectorXd::Zero(tot_dofs);
        f_S.segment(0, velocity_dofs_data.NumberDOFs) = f_x;
        f_S.segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs) = f_y;
      }

      auto initial_condition_function = [](const double &x, const double &y, const double &z) { return 0.0; };

      Eigen::VectorXd u_k = Eigen::VectorXd::Zero(tot_dofs);
      Eigen::VectorXd dp_strong = Eigen::VectorXd::Zero(pressure_dofs_data.NumberStrongs);
      Eigen::VectorXd du_x_strong = Eigen::VectorXd::Zero(velocity_dofs_data.NumberStrongs);
      Eigen::VectorXd du_y_strong = Eigen::VectorXd::Zero(velocity_dofs_data.NumberStrongs);
      Eigen::VectorXd u_x_numeric = u_k.segment(0, velocity_dofs_data.NumberDOFs);
      Eigen::VectorXd u_y_numeric = u_k.segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs);
      Eigen::VectorXd p_numeric = u_k.segment(2 * velocity_dofs_data.NumberDOFs, pressure_dofs_data.NumberDOFs);

      while (num_iteration < max_iterations && residual_norm > newton_tol * solution_norm)
      {
        const auto C = PDETools::Assembler_Utilities::PCC_2D::assemble_NS_operators(geometry_utilities,
                                                                                    mesh,
                                                                                    mesh_geometric_data,
                                                                                    velocity_dofs_data,
                                                                                    velocity_reference_element_data,
                                                                                    u_x_numeric,
                                                                                    u_y_numeric,
                                                                                    u_x_strong,
                                                                                    u_y_strong);

        Eigen::VectorXd dp_numeric;
        Eigen::VectorXd du_x_numeric;
        Eigen::VectorXd du_y_numeric;

        {
          const auto J_C = PDETools::Assembler_Utilities::PCC_2D::to_SparseMatrix(
                             PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(C.convective_operator.operator_dofs,
                                                                                         {tot_dofs, tot_dofs},
                                                                                         {0, 0}));

          Eigen::VectorXd du;
          {
            const auto J =
                PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(J_S + J_C);
            Eigen::VectorXd f_C = Eigen::VectorXd::Zero(tot_dofs);
            f_C.segment(0, 2 * velocity_dofs_data.NumberDOFs) = C.convective_rhs;

            auto f = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(f_S - f_C);
            const auto uk = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(u_k);
            f.SubtractionMultiplication(J_Stokes, uk);

            Gedim::Eigen_Array<> du_array;
            du_array.SetSize(tot_dofs);

            Gedim::Eigen_LUSolver solver;
            solver.Initialize(J);
            solver.Solve(f, du_array);
            du = PDETools::Assembler_Utilities::PCC_2D::to_VectorXd(du_array);
          }

          u_k = u_k + du;

          du_x_numeric = du.segment(0, velocity_dofs_data.NumberDOFs);
          du_y_numeric = du.segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs);
          dp_numeric = du.segment(2 * velocity_dofs_data.NumberDOFs, pressure_dofs_data.NumberDOFs);

          u_x_numeric = u_k.segment(0, velocity_dofs_data.NumberDOFs);
          u_y_numeric = u_k.segment(velocity_dofs_data.NumberDOFs, velocity_dofs_data.NumberDOFs);
          p_numeric = u_k.segment(2 * velocity_dofs_data.NumberDOFs, pressure_dofs_data.NumberDOFs);
        }

        const auto dp_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        pressure_dofs_data,
                                                                                        pressure_reference_element_data,
                                                                                        dp_numeric,
                                                                                        dp_strong);

        const auto du_x_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          velocity_dofs_data,
                                                                                          velocity_reference_element_data,
                                                                                          du_x_numeric,
                                                                                          du_x_strong);
        const auto du_y_error_L2 = PDETools::Assembler_Utilities::PCC_2D::compute_error_L2(geometry_utilities,
                                                                                          mesh,
                                                                                          mesh_geometric_data,
                                                                                          velocity_dofs_data,
                                                                                          velocity_reference_element_data,
                                                                                          du_y_numeric,
                                                                                          du_y_strong);

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

        solution_norm = std::sqrt(u_x_error_L2.numeric_norm_L2 *
                                  u_x_error_L2.numeric_norm_L2 +
                                  u_y_error_L2.numeric_norm_L2 *
                                  u_y_error_L2.numeric_norm_L2 +
                                  p_error_L2.numeric_norm_L2 *
                                  p_error_L2.numeric_norm_L2);
        residual_norm = std::sqrt(du_x_error_L2.numeric_norm_L2 *
                                  du_x_error_L2.numeric_norm_L2 +
                                  du_y_error_L2.numeric_norm_L2 *
                                  du_y_error_L2.numeric_norm_L2 +
                                  dp_error_L2.numeric_norm_L2 *
                                  dp_error_L2.numeric_norm_L2);
        num_iteration++;

        std::cout.precision(2);
        std::cout<< std::scientific<< "u_x_errorL2: "<< (u_x_error_L2.error_L2 / u_x_error_L2.numeric_norm_L2)<< " ";
        std::cout<< std::scientific<< "u_y_errorL2: "<< (u_y_error_L2.error_L2 / u_y_error_L2.numeric_norm_L2)<< " ";
        std::cout<< std::scientific<< "p_errorL2: "<< (p_error_L2.error_L2 / p_error_L2.numeric_norm_L2)<< std::endl;
        std::cout<< std::scientific<< "res: "<< (residual_norm / solution_norm)<< " / "<< newton_tol<< " ";
        std::cout<< " it: "<< num_iteration<< " / "<< max_iterations<< std::endl;
      }

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

      const auto u_x_on_quadrature =
          PDETools::Assembler_Utilities::PCC_2D::evaluate_solution_on_quadrature_points(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        velocity_dofs_data,
                                                                                        velocity_reference_element_data,
                                                                                        u_x_numeric,
                                                                                        u_x_strong,
                                                                                        velocity_x_exact_solution_function);
      const auto u_y_on_quadrature =
          PDETools::Assembler_Utilities::PCC_2D::evaluate_solution_on_quadrature_points(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        velocity_dofs_data,
                                                                                        velocity_reference_element_data,
                                                                                        u_y_numeric,
                                                                                        u_y_strong,
                                                                                        velocity_y_exact_solution_function);
      const auto p_on_quadrature =
          PDETools::Assembler_Utilities::PCC_2D::evaluate_solution_on_quadrature_points(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        pressure_dofs_data,
                                                                                        pressure_reference_element_data,
                                                                                        p_numeric,
                                                                                        p_strong,
                                                                                        pressure_exact_solution_function);

      const auto u_x_on_cell0Ds =
          PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh, velocity_dofs_data, u_x_numeric, u_x_strong, velocity_x_exact_solution_function);
      const auto u_y_on_cell0Ds =
          PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh, velocity_dofs_data, u_y_numeric, u_y_strong, velocity_y_exact_solution_function);
      const auto p_on_cell0Ds =
          PDETools::Assembler_Utilities::PCC_2D::extract_solution_on_cell0Ds(mesh, pressure_dofs_data, p_numeric, p_strong, pressure_exact_solution_function);

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"u_x_numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(u_x_on_cell0Ds.numeric_solution.size()),
                               u_x_on_cell0Ds.numeric_solution.data()},
                              {"u_x_exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(u_x_on_cell0Ds.exact_solution.size()),
                               u_x_on_cell0Ds.exact_solution.data()},
                              {"u_y_numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(u_y_on_cell0Ds.numeric_solution.size()),
                               u_y_on_cell0Ds.numeric_solution.data()},
                              {"u_y_exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(u_y_on_cell0Ds.exact_solution.size()),
                               u_y_on_cell0Ds.exact_solution.data()},
                              {"p_numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(p_on_cell0Ds.numeric_solution.size()),
                               p_on_cell0Ds.numeric_solution.data()},
                              {"p_exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(p_on_cell0Ds.exact_solution.size()),
                               p_on_cell0Ds.exact_solution.data()},
                              {"u_x_error_L2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(u_x_error_L2.cell2Ds_error_L2.size()),
                               u_x_error_L2.cell2Ds_error_L2.data()},
                              {"u_y_error_L2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(u_y_error_L2.cell2Ds_error_L2.size()),
                               u_y_error_L2.cell2Ds_error_L2.data()},
                              {"p_error_L2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(p_error_L2.cell2Ds_error_L2.size()),
                               p_error_L2.cell2Ds_error_L2.data()}});
        exporter.Export(exportFolder + "/solution.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPoints(u_x_on_quadrature.quadrature_points,
                           {{"u_x_numeric",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(u_x_on_quadrature.numeric_solution.size()),
                             u_x_on_quadrature.numeric_solution.data()},
                            {"u_y_numeric",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(u_y_on_quadrature.numeric_solution.size()),
                             u_y_on_quadrature.numeric_solution.data()},
                            {"u_x_exact",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(u_x_on_quadrature.exact_solution.size()),
                             u_x_on_quadrature.exact_solution.data()},
                            {"u_y_exact",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(u_y_on_quadrature.exact_solution.size()),
                             u_y_on_quadrature.exact_solution.data()}});

        exporter.Export(exportFolder + "/u_on_quadrature.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPoints(p_on_quadrature.quadrature_points,
                           {

                             {"p_numeric",
                              Gedim::VTPProperty::Formats::Points,
                              static_cast<unsigned int>(p_on_quadrature.numeric_solution.size()),
                              p_on_quadrature.numeric_solution.data()},
                             {"p_exact",
                              Gedim::VTPProperty::Formats::Points,
                              static_cast<unsigned int>(p_on_quadrature.exact_solution.size()),
                              p_on_quadrature.exact_solution.data()}});

        exporter.Export(exportFolder + "/p_on_quadrature.vtu");
      }

      ASSERT_TRUE(u_x_error_L2.error_L2 < 1.0e-1 * u_x_error_L2.numeric_norm_L2);
      ASSERT_TRUE(u_y_error_L2.error_L2 < 1.0e-1 * u_y_error_L2.numeric_norm_L2);
      ASSERT_TRUE(p_error_L2.error_L2 < 1.0e-1 * p_error_L2.numeric_norm_L2);
    }

  } // namespace UnitTesting
} // namespace Polydim

#endif
