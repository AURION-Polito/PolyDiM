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

#ifndef __TEST_assembler_PCC_2D_H
#define __TEST_assembler_PCC_2D_H

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

    TEST(TEST_assembler_PCC_2D, TEST_assembler_PCC_2D_forcing_term)
    {
      Gedim::GeometryUtilitiesConfig geometry_utilities_config;
      geometry_utilities_config.Tolerance1D = 1.0e-8;
      geometry_utilities_config.Tolerance2D = 1.0e-12;
      Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

      const std::string exportFolder = "./Export/TEST_assembler_PCC_2D/TEST_assembler_PCC_2D_forcing_term";
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
                                                                  0.1,
                                                                  mesh);

      const double a = 2.0;
      const std::array<double, 3> b = { 0.0, 0.0, 0.0 };
      const double c = 3.0;
      const unsigned int method_order = 3;
      const auto reference_element_data =
          Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC,
                                                                       method_order);

      const auto mesh_geometric_data = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(
                                         geometry_utilities,
                                         mesh_utilities,
                                         mesh,
                                         Polydim::PDETools::LocalSpace_PCC_2D::MeshGeometricDataConfigiguration(reference_element_data));

      Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);
      Polydim::PDETools::DOFs::DOFsManager dofManager;

      std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info = {
        {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
        {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
        {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
        {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
        {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

      const auto mesh_dofs_info = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
      const auto dofs_data = dofManager.CreateDOFs_2D(mesh_dofs_info, mesh_connectivity_data);

      auto exact_solution_function = [&method_order, &a, &b, &c](const double &x, const double &y, const double &z) {
        const double polynomial = x + y + 0.5;

        double result = 1.0;
        for (int i = 0; i < method_order; ++i)
          result *= polynomial;

        return result;
      };

      auto exact_gradient_solution_function = [&method_order, &a, &b, &c](const double &x, const double &y, const double &z) {
        const double polynomial = x + y + 0.5;

        double result = method_order;
        const int max_order = method_order - 1;
        for (int i = 0; i < max_order; ++i)
          result *= polynomial;

        return std::array<double, 3>({ result, result, 0.0 });
      };

      auto exact_laplacian_solution_function = [&method_order, &a, &b, &c](const double &x, const double &y, const double &z) {
        double result = 2.0 * method_order * (method_order - 1);
        const double polynomial = x + y + 0.5;

        const int max_order = method_order - 2;
        for (int i = 0; i < max_order; ++i)
          result *= polynomial;

        return result;
      };

      auto source_term_function = [&method_order, &a, &b, &c, &exact_solution_function, &exact_gradient_solution_function, &exact_laplacian_solution_function](const double &x, const double &y, const double &z, const Eigen::VectorXd &check) {
        const auto u_lap = exact_laplacian_solution_function(x, y, z);
        const auto u_grad = exact_gradient_solution_function(x, y, z);
        const auto u = exact_solution_function(x, y, z);

        return - a * u_lap +
            (b.at(0) * u_grad.at(0) +
             b.at(1) * u_grad.at(1) +
             b.at(2) * u_grad.at(2)) +
            c * u;
      };

      const auto source_term = PDETools::Assembler_Utilities::PCC_2D::assembler_source_term(geometry_utilities,
                                                                                            mesh,
                                                                                            mesh_geometric_data,
                                                                                            dofs_data,
                                                                                            reference_element_data,
                                                                                            source_term_function);

      ASSERT_EQ(dofs_data.NumberDOFs, source_term.size());

      auto diffusion_term_function = [&a](const double &x, const double &y, const double &z) {
        return a;
      };

      auto advection_term_function = [&b](const double &x, const double &y, const double &z) {
        return b;
      };

      auto reaction_term_function = [&c](const double &x, const double &y, const double &z) {
        return c;
      };

      const auto elliptic_operator = PDETools::Assembler_Utilities::PCC_2D::assembler_elliptic_operator(geometry_utilities,
                                                                                                        mesh,
                                                                                                        mesh_geometric_data,
                                                                                                        dofs_data,
                                                                                                        reference_element_data,
                                                                                                        diffusion_term_function,
                                                                                                        reaction_term_function);

      ASSERT_EQ(dofs_data.NumberDOFs, elliptic_operator.A.size.at(0));
      ASSERT_EQ(dofs_data.NumberDOFs, elliptic_operator.A.size.at(1));
      ASSERT_EQ(dofs_data.NumberDOFs, elliptic_operator.A_Strong.size.at(0));
      ASSERT_EQ(dofs_data.NumberStrongs, elliptic_operator.A_Strong.size.at(1));

      auto strong_solution_function =
          [&exact_solution_function](const unsigned int marker, const double &x, const double &y, const double &z) {
        if (marker != 1)
          throw std::runtime_error("marker not managed");

        return exact_solution_function(x, y, z);
      };

      const auto strong_solution = PDETools::Assembler_Utilities::PCC_2D::assembler_strong_solution(geometry_utilities,
                                                                                                    mesh,
                                                                                                    mesh_geometric_data,
                                                                                                    mesh_dofs_info,
                                                                                                    dofs_data,
                                                                                                    reference_element_data,
                                                                                                    strong_solution_function);

      ASSERT_EQ(dofs_data.NumberStrongs, strong_solution.size());

      auto weak_term_function = [&a, &exact_gradient_solution_function](const unsigned int marker, const double &x, const double &y, const double &z) {

        const auto u_grad = exact_gradient_solution_function(x, y, z);

        switch (marker)
        {
          case 2:
            return a * u_grad.at(0);
          case 4:
            return a * u_grad.at(1);
          default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
      };

      const auto weak_term = PDETools::Assembler_Utilities::PCC_2D::assembler_weak_term(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        mesh_dofs_info,
                                                                                        dofs_data,
                                                                                        reference_element_data,
                                                                                        weak_term_function);

      ASSERT_EQ(dofs_data.NumberDOFs, weak_term.size());

      const auto exact_solution = PDETools::Assembler_Utilities::PCC_2D::assembler_exact_solution(geometry_utilities,
                                                                                                  mesh,
                                                                                                  mesh_geometric_data,
                                                                                                  dofs_data,
                                                                                                  reference_element_data,
                                                                                                  exact_solution_function);

      ASSERT_EQ(dofs_data.NumberDOFs, exact_solution.exact_solution.size());
      ASSERT_EQ(dofs_data.NumberStrongs, exact_solution.exact_solution_strong.size());

      {
        const auto f = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(source_term);
        const auto w_t = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(weak_term);
        const auto u_D = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_Array(strong_solution);
        const auto A = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(elliptic_operator.A);
        const auto A_D = PDETools::Assembler_Utilities::PCC_2D::to_Eigen_SparseArray(elliptic_operator.A_Strong);

        auto rhs = f;
        rhs += w_t;
        rhs.SubtractionMultiplication(A_D, u_D);

        Gedim::Eigen_Array<> u;
        u.SetSize(dofs_data.NumberDOFs);

        Gedim::Eigen_LUSolver solver;
        solver.Initialize(A);
        solver.Solve(rhs, u);

        const auto numeric_solution = PDETools::Assembler_Utilities::PCC_2D::to_VectorXd(u);

        const auto u_cell0Ds = PDETools::Assembler_Utilities::PCC_2D::assembler_extract_cell0Ds(mesh,
                                                                                        dofs_data,
                                                                                        numeric_solution,
                                                                                        strong_solution,
                                                                                        exact_solution_function);


        const auto error_L2 = PDETools::Assembler_Utilities::PCC_2D::assembler_error_L2(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        dofs_data,
                                                                                        reference_element_data,
                                                                                        numeric_solution,
                                                                                        strong_solution,
                                                                                        exact_solution_function);

        const auto error_H1 = PDETools::Assembler_Utilities::PCC_2D::assembler_error_H1(geometry_utilities,
                                                                                        mesh,
                                                                                        mesh_geometric_data,
                                                                                        dofs_data,
                                                                                        reference_element_data,
                                                                                        numeric_solution,
                                                                                        strong_solution,
                                                                                        exact_gradient_solution_function);

        {
          Gedim::VTKUtilities exporter;

          exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                               mesh.Cell2DsVertices(),
                               {
                                 {
                                  "numeric_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(u_cell0Ds.cell0Ds_numeric.size()),
                                   u_cell0Ds.cell0Ds_numeric.data()
                                 },
                                 {
                                  "exact_solution",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(u_cell0Ds.cell0Ds_exact.size()),
                                   u_cell0Ds.cell0Ds_exact.data()
                                 },
                                 {
                                  "error_L2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(error_L2.cell2Ds_error_L2.size()),
                                   error_L2.cell2Ds_error_L2.data()
                                 },
                                 {
                                  "error_H1",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(error_H1.cell2Ds_error_H1.size()),
                                   error_H1.cell2Ds_error_H1.data()
                                 }
                               });
          exporter.Export(exportFolder + "/solution.vtu");
        }

        std::cout.precision(2);
        //std::cout<< std::scientific<< "A: "<< A<< std::endl;
        //std::cout<< std::scientific<< "A_D: "<< A_D<< std::endl;
        std::cout << std::scientific << "f: " << f << std::endl;
        std::cout << std::scientific << "w_t: " << w_t << std::endl;
        std::cout << std::scientific << "r: " << rhs << std::endl;
        std::cout << std::scientific << "u: " << u << std::endl;
        std::cout << std::scientific << "u_ex: " << exact_solution.exact_solution.transpose() << std::endl;
        std::cout << std::scientific << "u_D: " << u_D << std::endl;
        std::cout << std::scientific << "u_ex_D: " << exact_solution.exact_solution_strong.transpose() << std::endl;
        std::cout << std::scientific << "err_L2: " << (error_L2.error_L2 / error_L2.exact_norm_L2) << std::endl;
        std::cout << std::scientific << "err_H1: " << (error_H1.error_H1 / error_H1.exact_norm_H1) << std::endl;

        ASSERT_TRUE((strong_solution - exact_solution.exact_solution_strong).norm() <
                    1.0e-13 * exact_solution.exact_solution_strong.norm());
        ASSERT_TRUE((numeric_solution - exact_solution.exact_solution).norm() < 1.0e-13 * exact_solution.exact_solution.norm());
        ASSERT_TRUE(error_L2.error_L2 < 1.0e-13 * error_L2.exact_norm_L2);
        ASSERT_TRUE(error_H1.error_H1 < 1.0e-13 * error_H1.exact_norm_H1);
      }
    }

  } // namespace UnitTesting
} // namespace Polydim

#endif
