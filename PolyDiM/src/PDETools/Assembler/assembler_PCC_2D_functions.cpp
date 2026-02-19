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

#include "assembler_PCC_2D_functions.hpp"
#include "EllipticEquation.hpp"
#include "assembler_PCC_2D_functions_utilities.hpp"

namespace Polydim
{
  namespace PDETools
  {
    namespace Assembler_Utilities
    {
      namespace PCC_2D
      {
        // ***************************************************************************
        Variational_Operator convert_operator(const Gedim::Eigen_SparseArray<> &elliptic_matrix, const Gedim::Eigen_SparseArray<> &elliptic_strong_matrix)
        {
          Variational_Operator result;

          result.A = to_Sparse_Matrix_Data(elliptic_matrix);
          result.A_Strong = to_Sparse_Matrix_Data(elliptic_strong_matrix);

          return result;
        }
        // ***************************************************************************
        Eigen::VectorXd function_evaluation(const Eigen::MatrixXd &points,
                                            const std::function<double(const double &, const double &, const double &, const Eigen::VectorXd &)> f)
        {
          Eigen::VectorXd function_values(points.cols());

          for (int i = 0; i < points.cols(); ++i)
          {
            function_values[i] = f(points(0, i), points(1, i), points(2, i), function_values);
          }

          return function_values;
        }
        // ***************************************************************************
        Eigen::VectorXd function_evaluation(const unsigned int marker,
                                            const Eigen::MatrixXd &points,
                                            const std::function<double(const unsigned int, const double &, const double &, const double &)> f)
        {
          Eigen::VectorXd function_values(points.cols());

          for (int i = 0; i < points.cols(); ++i)
          {
            function_values[i] = f(marker, points(0, i), points(1, i), points(2, i));
          }

          return function_values;
        }
        // ***************************************************************************
        Eigen::VectorXd function_evaluation(const Eigen::MatrixXd &points,
                                            const std::function<double(const double &, const double &, const double &)> f)
        {
          Eigen::VectorXd function_values(points.cols());

          for (int i = 0; i < points.cols(); ++i)
          {
            function_values[i] = f(points(0, i), points(1, i), points(2, i));
          }

          return function_values;
        }
        // ***************************************************************************
        std::array<Eigen::VectorXd, 3> function_evaluation(const Eigen::MatrixXd &points,
                                                           const std::function<std::array<double, 3>(const double &, const double &, const double &)> f)
        {
          std::array<Eigen::VectorXd, 3> function_values;
          function_values.at(0).resize(points.cols());
          function_values.at(1).resize(points.cols());
          function_values.at(2).resize(points.cols());

          for (int i = 0; i < points.cols(); ++i)
          {
            const auto result_f = f(points(0, i), points(1, i), points(2, i));

            function_values.at(0)[i] = result_f.at(0);
            function_values.at(1)[i] = result_f.at(1);
            function_values.at(2)[i] = result_f.at(2);
          }

          return function_values;
        }
        // ***************************************************************************
        Eigen::VectorXd assembler_source_term(const Gedim::GeometryUtilities &geometry_utilities,
                                              const Gedim::MeshMatricesDAO &mesh,
                                              const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                              const DOFs::DOFsManager::DOFsData &test_dofs_data,
                                              const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                              const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &test_reference_element_data,
                                              const std::function<double(const double &, const double &, const double &, const Eigen::VectorXd &)> source_term_function)
        {
          Gedim::Eigen_Array<> forcing_term;

          forcing_term.SetSize(test_dofs_data.NumberDOFs);

          Polydim::PDETools::Equations::EllipticEquation equation;

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(test_dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                      geometry_utilities.Tolerance2D(),
                                                                                                      mesh_geometric_data,
                                                                                                      c,
                                                                                                      trial_reference_element_data);

            const auto test_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                      geometry_utilities.Tolerance2D(),
                                                                                                      mesh_geometric_data,
                                                                                                      c,
                                                                                                      test_reference_element_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(trial_reference_element_data, trial_local_space_data);

            const auto test_basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(test_reference_element_data,
                                                                           test_local_space_data,
                                                                           cell2D_internal_quadrature.Points);

            const Eigen::VectorXd source_term_values = function_evaluation(cell2D_internal_quadrature.Points, source_term_function);

            const Eigen::VectorXd local_rhs =
                equation.ComputeCellForcingTerm(source_term_values, test_basis_functions_values, cell2D_internal_quadrature.Weights);

            assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(test_reference_element_data, test_local_space_data) ==
                   test_dofs_data.CellsGlobalDOFs[2].at(c).size());

            Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                              local_matrix_to_global_matrix_dofs_data,
                                                                                              local_rhs,
                                                                                              forcing_term);
          }

          forcing_term.Create();

          return static_cast<Eigen::VectorXd &>(forcing_term);
        }
        // ***************************************************************************
        Variational_Operator assemble_elliptic_operator(const Gedim::GeometryUtilities &geometry_utilities,
                                                         const Gedim::MeshMatricesDAO &mesh,
                                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                         const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                         const DOFs::DOFsManager::DOFsData &test_dofs_data,
                                                         const LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                                         const LocalSpace_PCC_2D::ReferenceElement_Data &test_reference_element_data,
                                                         const std::function<double(const double&, const double&, const double&)> diffusion_term_function,
                                                         const std::function<std::array<double, 3>(const double&, const double&, const double&)> advection_term_function,
                                                         const std::function<double(const double &, const double &, const double &)> reaction_term_function)
        {
          Gedim::Eigen_SparseArray<> elliptic_matrix;
          Gedim::Eigen_SparseArray<> elliptic_strong_matrix;

          elliptic_matrix.SetSize(test_dofs_data.NumberDOFs,
                                  trial_dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::None);
          elliptic_strong_matrix.SetSize(test_dofs_data.NumberDOFs,
                                         trial_dofs_data.NumberStrongs);

          Polydim::PDETools::Equations::EllipticEquation equation;

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data trial_local_matrix_to_global_matrix_dofs_data =
          {{std::cref(trial_dofs_data)}, {0}, {0}, {0}};
          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data test_local_matrix_to_global_matrix_dofs_data =
          {{std::cref(test_dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                       geometry_utilities.Tolerance2D(),
                                                                                                       mesh_geometric_data,
                                                                                                       c,
                                                                                                       trial_reference_element_data);


            const auto test_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                      geometry_utilities.Tolerance2D(),
                                                                                                      mesh_geometric_data,
                                                                                                      c,
                                                                                                      test_reference_element_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(trial_reference_element_data, trial_local_space_data);


            const auto trial_basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(trial_reference_element_data, trial_local_space_data);
            const auto trial_basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(trial_reference_element_data, trial_local_space_data);

            const auto test_basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(test_reference_element_data,
                                                                           test_local_space_data,
                                                                           cell2D_internal_quadrature.Points);
            const auto test_basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(test_reference_element_data,
                                                                                     test_local_space_data,
                                                                                     cell2D_internal_quadrature.Points);

            const auto test_local_space_size = Polydim::PDETools::LocalSpace_PCC_2D::Size(test_reference_element_data, test_local_space_data);
            const auto trial_local_space_size = Polydim::PDETools::LocalSpace_PCC_2D::Size(trial_reference_element_data, trial_local_space_data);

            Eigen::MatrixXd local_A = Eigen::MatrixXd::Zero(test_local_space_size,
                                                            trial_local_space_size);

            if (diffusion_term_function)
            {
              const auto diffusion_term_values = function_evaluation(cell2D_internal_quadrature.Points, diffusion_term_function);

              local_A += equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                             trial_basis_functions_derivative_values,
                                                             test_basis_functions_derivative_values,
                                                             cell2D_internal_quadrature.Weights);

              if (test_local_space_size == trial_local_space_size)
              {
                const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
                local_A +=
                    k_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(test_reference_element_data, test_local_space_data);
              }
            }

            if (advection_term_function)
            {
              const auto advection_term_values = function_evaluation(cell2D_internal_quadrature.Points,
                                                                     advection_term_function);

              local_A += equation.ComputeCellAdvectionMatrix(advection_term_values,
                                                             test_basis_functions_values,
                                                             trial_basis_functions_derivative_values,
                                                             cell2D_internal_quadrature.Weights);

            }

            if (reaction_term_function)
            {
              const auto reaction_term_values = function_evaluation(cell2D_internal_quadrature.Points, reaction_term_function);

              local_A += equation.ComputeCellReactionMatrix(reaction_term_values,
                                                            trial_basis_functions_values,
                                                            test_basis_functions_values,
                                                            cell2D_internal_quadrature.Weights);
            }

            assert(trial_local_space_size ==
                   trial_dofs_data.CellsGlobalDOFs[2].at(c).size());
            assert(test_local_space_size ==
                   test_dofs_data.CellsGlobalDOFs[2].at(c).size());

            Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                              test_local_matrix_to_global_matrix_dofs_data,
                                                                                              trial_local_matrix_to_global_matrix_dofs_data,
                                                                                              local_A,
                                                                                              elliptic_matrix,
                                                                                              elliptic_strong_matrix);
          }

          elliptic_matrix.Create();
          elliptic_strong_matrix.Create();

          return convert_operator(elliptic_matrix, elliptic_strong_matrix);
        }
        // ***************************************************************************
        Eigen::VectorXd assemble_strong_solution(const Gedim::GeometryUtilities &geometry_utilities,
                                                  const Gedim::MeshMatricesDAO &mesh,
                                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                  const DOFs::DOFsManager::MeshDOFsInfo &trial_mesh_dofs_info,
                                                  const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                  const LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                                  const std::function<double(const unsigned int, const double &, const double &, const double &)> strong_solution_function)
        {
          Gedim::Eigen_Array<> strong_solution;

          strong_solution.SetSize(trial_dofs_data.NumberStrongs);

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(trial_dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); ++v)
            {
              const unsigned int cell0D_index = mesh.Cell2DVertex(c, v);
              const auto &boundary_info = trial_mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                continue;

              const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

              const auto strong_boundary_values = function_evaluation(boundary_info.Marker, coordinates, strong_solution_function);

              const auto local_dofs = trial_dofs_data.CellsDOFs.at(0).at(cell0D_index);

              assert(local_dofs.size() == static_cast<unsigned int>(strong_boundary_values.size()));

              for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
              {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    strong_solution.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
                  }
                    break;
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    continue;
                  default:
                    throw std::runtime_error("Unknown DOF Type");
                }
              }
            }

            if (trial_reference_element_data.Order > 1)
            {
              const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                         geometry_utilities.Tolerance2D(),
                                                                                                         mesh_geometric_data,
                                                                                                         c,
                                                                                                         trial_reference_element_data);

              for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
              {
                const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

                const auto &boundary_info = trial_mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
                const auto local_dofs = trial_dofs_data.CellsDOFs.at(1).at(cell1D_index);

                if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
                    local_dofs.size() == 0)
                  continue;

                const auto edge_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(trial_reference_element_data, trial_local_space_data, ed);

                const auto strong_boundary_values =
                    function_evaluation(boundary_info.Marker, edge_dofs_coordinates, strong_solution_function);

                assert(local_dofs.size() == static_cast<unsigned int>(strong_boundary_values.size()));

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                  const auto &local_dof_i = local_dofs.at(loc_i);

                  switch (local_dof_i.Type)
                  {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                      strong_solution.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
                    }
                      break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                      continue;
                    default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
                }
              }

            }
          }

          strong_solution.Create();

          return static_cast<Eigen::VectorXd &>(strong_solution);
        }
        // ***************************************************************************
        Exact_Solution_Data assemble_exact_solution(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                                     const std::function<double(const double &, const double &, const double &)> exact_solution_function)
        {
          Gedim::Eigen_Array<> exact_solution;
          Gedim::Eigen_Array<> exact_solution_strong;

          exact_solution.SetSize(trial_dofs_data.NumberDOFs);
          exact_solution_strong.SetSize(trial_dofs_data.NumberStrongs);

          // Assemble equation elements
          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
          {
            // DOFs: vertices
            const Eigen::MatrixXd coordinates = mesh.Cell2DVerticesCoordinates(c);
            const Eigen::VectorXd dofs_vertices = function_evaluation(coordinates, exact_solution_function);

            // Assemble local numerical solution
            unsigned int count = 0;
            for (unsigned int p = 0; p < mesh.Cell2DNumberVertices(c); p++)
            {
              const unsigned int cell0D_index = mesh.Cell2DVertex(c, p);

              const auto local_dofs = trial_dofs_data.CellsDOFs.at(0).at(cell0D_index);
              for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
              {
                const auto &local_dof_i = local_dofs.at(loc_i);
                const int global_i = local_dof_i.Global_Index;

                switch (local_dof_i.Type)
                {
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    exact_solution_strong.SetValue(global_i, dofs_vertices(count++));
                  }
                    break;
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    exact_solution.SetValue(global_i, dofs_vertices(count++));
                  }
                    break;
                  default:
                    throw std::runtime_error("Unknown DOF Type");
                }
              }
            }

            // Assemble strong boundary condition on Cell1Ds
            if (trial_reference_element_data.Order > 1)
            {
              const auto trial_local_space_data =
                  Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                         geometry_utilities.Tolerance2D(),
                                                                         mesh_geometric_data,
                                                                         c,
                                                                         trial_reference_element_data);

              // Assemble strong boundary condition on Cell1Ds
              for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
              {
                const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

                const auto local_dofs = trial_dofs_data.CellsDOFs.at(1).at(cell1D_index);

                const auto edge_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(trial_reference_element_data, trial_local_space_data, ed);

                const Eigen::VectorXd dofs_edge = function_evaluation(edge_dofs_coordinates, exact_solution_function);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                  const auto &local_dof_i = local_dofs.at(loc_i);
                  const int global_i = local_dof_i.Global_Index;

                  switch (local_dof_i.Type)
                  {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                      exact_solution_strong.SetValue(global_i, dofs_edge(loc_i));
                    }
                      break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                      exact_solution.SetValue(global_i, dofs_edge(loc_i));
                    }
                      break;
                    default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
                }
              }

              const auto local_dofs = trial_dofs_data.CellsDOFs.at(2).at(c);

              if (local_dofs.size())
              {
                const auto internal_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofsCoordinates(trial_reference_element_data, trial_local_space_data);

                const Eigen::VectorXd exact_values_at_dofs =
                    function_evaluation(internal_dofs_coordinates.Points, exact_solution_function);

                const Eigen::VectorXd dofs_internal =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofs(trial_reference_element_data,
                                                                       trial_local_space_data,
                                                                       exact_values_at_dofs,
                                                                       internal_dofs_coordinates);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                  const auto &local_dof_i = local_dofs.at(loc_i);
                  const int global_i = local_dof_i.Global_Index;

                  switch (local_dof_i.Type)
                  {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                      exact_solution_strong.SetValue(global_i, dofs_internal(loc_i));
                    }
                      break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                      exact_solution.SetValue(global_i, dofs_internal(loc_i));
                    }
                      break;
                    default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
                }
              }
            }
          }

          exact_solution.Create();
          exact_solution_strong.Create();

          return {static_cast<Eigen::VectorXd &>(exact_solution), static_cast<Eigen::VectorXd &>(exact_solution_strong)};
        }
        // ***************************************************************************
        Eigen::VectorXd assemble_weak_term(const Gedim::GeometryUtilities &geometry_utilities,
                                            const Gedim::MeshMatricesDAO &mesh,
                                            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                            const DOFs::DOFsManager::MeshDOFsInfo &trial_mesh_dofs_info,
                                            const DOFs::DOFsManager::DOFsData &test_dofs_data,
                                            const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                            const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &test_reference_element_data,
                                            const std::function<double(const unsigned int, const double &, const double &, const double &)> weak_term_function)
        {
          Gedim::Eigen_Array<> weak_term;

          weak_term.SetSize(test_dofs_data.NumberDOFs);

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;


            const auto test_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 test_reference_element_data);

            const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(c).cols();

            for (unsigned int ed = 0; ed < numVertices; ed++)
            {
              const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

              const auto &boundary_info = trial_mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
                continue;

              // compute vem values
              const auto weakReferenceSegment =
                  Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * trial_reference_element_data.Order);

              const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

              // map edge internal quadrature points
              const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell2DsEdgeDirections.at(c)[ed]
                                                 ? mesh_geometric_data.Cell2DsVertices.at(c).col(ed)
                                                 : mesh_geometric_data.Cell2DsVertices.at(c).col((ed + 1) % numVertices);

              const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell2DsEdgeTangents.at(c).col(ed);
              const double direction = mesh_geometric_data.Cell2DsEdgeDirections.at(c)[ed] ? 1.0 : -1.0;

              const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
              Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
              for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
                weakQuadraturePoints.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;

              const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(c)[ed]);
              const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

              const Eigen::VectorXd neumannValues = function_evaluation(boundary_info.Marker, weakQuadraturePoints, weak_term_function);
              const auto weak_basis_function_values =
                  Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValuesOnEdge(ed, test_reference_element_data, test_local_space_data, pointsCurvilinearCoordinates);

              // compute values of Neumann condition
              const Eigen::VectorXd neumannContributions =
                  weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues;

              for (unsigned int p = 0; p < 2; ++p)
              {
                const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index, p);

                const auto local_dofs = test_dofs_data.CellsDOFs.at(0).at(cell0D_index);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                  const auto &local_dof_i = local_dofs.at(loc_i);

                  switch (local_dof_i.Type)
                  {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                      continue;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                      weak_term.AddValue(local_dof_i.Global_Index, neumannContributions[p]);
                    }
                      break;
                    default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
                }
              }

              const auto local_dofs = test_dofs_data.CellsDOFs.at(1).at(cell1D_index);
              for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
              {
                const auto &local_dof_i = local_dofs.at(loc_i);

                const unsigned int localIndex = loc_i;

                switch (local_dof_i.Type)
                {
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    continue;
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    weak_term.AddValue(local_dof_i.Global_Index, neumannContributions[localIndex + 2]);
                  }
                    break;
                  default:
                    throw std::runtime_error("Unknown DOF Type");
                }
              }
            }
          }

          weak_term.Create();

          return static_cast<Eigen::VectorXd &>(weak_term);
        }
        // ***************************************************************************
        Post_Process_Data_Cell0Ds extract_solution_on_cell0Ds(const Gedim::MeshMatricesDAO &mesh,
                                                            const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                            const Eigen::VectorXd &numerical_solution,
                                                            const Eigen::VectorXd &numerical_solution_strong,
                                                            const std::function<double(const double &, const double &, const double &)> exact_solution_function,
                                                              const std::function<std::array<double, 3>(const double &, const double &, const double &)> exact_gradient_solution_function)
        {
          Post_Process_Data_Cell0Ds result;

          const auto num_solution = to_Eigen_Array(numerical_solution);
          const auto num_solution_strong = to_Eigen_Array(numerical_solution_strong);

          result.numeric_solution.resize(mesh.Cell0DTotalNumber());
          result.exact_solution.resize(mesh.Cell0DTotalNumber());
          result.exact_gradient_solution.at(0).resize(mesh.Cell0DTotalNumber());
          result.exact_gradient_solution.at(1).resize(mesh.Cell0DTotalNumber());
          result.exact_gradient_solution.at(2).resize(mesh.Cell0DTotalNumber());

          for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
          {
            if (!mesh.Cell0DIsActive(p))
              continue;

            result.exact_solution[p] = function_evaluation(mesh.Cell0DCoordinates(p), exact_solution_function)[0];
            const auto grad_solution = function_evaluation(mesh.Cell0DCoordinates(p), exact_gradient_solution_function);

            result.exact_gradient_solution.at(0)[p] = grad_solution.at(0)[0];
            result.exact_gradient_solution.at(1)[p] = grad_solution.at(1)[0];
            result.exact_gradient_solution.at(2)[p] = grad_solution.at(2)[0];

            const auto local_dofs = trial_dofs_data.CellsDOFs.at(0).at(p);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
              const auto &local_dof_i = local_dofs.at(loc_i);

              switch (local_dof_i.Type)
              {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                  result.numeric_solution[p] = numerical_solution_strong[local_dof_i.Global_Index];
                  break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                  result.numeric_solution[p] = numerical_solution[local_dof_i.Global_Index];
                  break;
                default:
                  throw std::runtime_error("Unknown DOF Type");
              }
            }
          }

          return result;
        }
        // ***************************************************************************
        Post_Process_Data_ErrorL2 compute_error_L2(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                                     const Eigen::VectorXd &numerical_solution,
                                                     const Eigen::VectorXd &numerical_solution_strong,
                                                     const std::function<double(const double &, const double &, const double &)> exact_solution_function)
        {
          Post_Process_Data_ErrorL2 result;

          const auto num_solution = to_Eigen_Array(numerical_solution);
          const auto num_solution_strong = to_Eigen_Array(numerical_solution_strong);

          result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
          result.cell2Ds_exact_norm_L2.setZero(mesh.Cell2DTotalNumber());
          result.cell2Ds_numeric_norm_L2.setZero(mesh.Cell2DTotalNumber());
          result.error_L2 = 0.0;
          result.numeric_norm_L2 = 0.0;
          result.exact_norm_L2 = 0.0;

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 trial_reference_element_data);

            const auto trial_basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(trial_reference_element_data,
                                                                           trial_local_space_data,
                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0k);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(trial_reference_element_data, trial_local_space_data);

            const auto exact_solution_values = function_evaluation(cell2D_internal_quadrature.Points, exact_solution_function);

            const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, trial_dofs_data);
            const Eigen::VectorXd dofs_values =
                PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                    trial_dofs_data,
                                                                                    local_count_dofs.num_total_dofs,
                                                                                    local_count_dofs.offsets_DOFs,
                                                                                    {0},
                                                                                    {0},
                                                                                    num_solution,
                                                                                    num_solution_strong);

            const Eigen::VectorXd local_error_L2 = (trial_basis_functions_values * dofs_values - exact_solution_values).array().square();
            const Eigen::VectorXd local_numeric_norm_L2 = (trial_basis_functions_values * dofs_values).array().square();
            const Eigen::VectorXd local_exact_norm_L2 = (exact_solution_values).array().square();

            result.cell2Ds_error_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2;
            result.cell2Ds_numeric_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_numeric_norm_L2;
            result.cell2Ds_exact_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_exact_norm_L2;
          }

          result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
          result.numeric_norm_L2 = std::sqrt(result.cell2Ds_numeric_norm_L2.sum());
          result.exact_norm_L2 = std::sqrt(result.cell2Ds_exact_norm_L2.sum());

          return result;
        }
        // ***************************************************************************
        Post_Process_Data_ErrorH1 compute_error_H1(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &trial_dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &trial_reference_element_data,
                                                     const Eigen::VectorXd &numerical_solution,
                                                     const Eigen::VectorXd &numerical_solution_strong,
                                                     const std::function<std::array<double, 3>(const double &, const double &, const double &)> exact_gradient_solution_function)
        {
          Post_Process_Data_ErrorH1 result;

          const auto num_solution = to_Eigen_Array(numerical_solution);
          const auto num_solution_strong = to_Eigen_Array(numerical_solution_strong);

          result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
          result.cell2Ds_exact_norm_H1.setZero(mesh.Cell2DTotalNumber());
          result.cell2Ds_numeric_norm_H1.setZero(mesh.Cell2DTotalNumber());
          result.error_H1 = 0.0;
          result.numeric_norm_H1 = 0.0;
          result.exact_norm_H1 = 0.0;

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 trial_reference_element_data);

            const auto trial_basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(trial_reference_element_data, trial_local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(trial_reference_element_data, trial_local_space_data);

            const auto exact_derivative_solution_values = function_evaluation(cell2D_internal_quadrature.Points, exact_gradient_solution_function);

            const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, trial_dofs_data);
            const Eigen::VectorXd dofs_values =
                PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                    trial_dofs_data,
                                                                                    local_count_dofs.num_total_dofs,
                                                                                    local_count_dofs.offsets_DOFs,
                                                                                    {0},
                                                                                    {0},
                                                                                    num_solution,
                                                                                    num_solution_strong);

            const Eigen::VectorXd local_error_H1 =
                (trial_basis_functions_derivative_values.at(0) * dofs_values -
                 exact_derivative_solution_values.at(0)).array().square() +
                (trial_basis_functions_derivative_values.at(1) * dofs_values -
                 exact_derivative_solution_values.at(1)).array().square();

            const Eigen::VectorXd local_numeric_norm_H1 =
                (trial_basis_functions_derivative_values.at(0) * dofs_values).array().square() +
                (trial_basis_functions_derivative_values.at(1) * dofs_values).array().square();

            const Eigen::VectorXd local_exact_norm_H1 =
                (exact_derivative_solution_values.at(0)).array().square() +
                (exact_derivative_solution_values.at(1)).array().square();

            result.cell2Ds_error_H1[c] = cell2D_internal_quadrature.Weights.transpose() *
                                         local_error_H1;
            result.cell2Ds_numeric_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_numeric_norm_H1;
            result.cell2Ds_exact_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_exact_norm_H1;
          }

          result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
          result.numeric_norm_H1 = std::sqrt(result.cell2Ds_numeric_norm_H1.sum());
          result.exact_norm_H1 = std::sqrt(result.cell2Ds_exact_norm_H1.sum());

          return result;
        }
// ***************************************************************************
        Evaluate_Solution_On_Quadrature_Points_Data evaluate_solution_on_quadrature_points(const Gedim::GeometryUtilities& geometry_utilities,
                                                                                           const Gedim::MeshMatricesDAO& mesh,
                                                                                           const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                                                           const DOFs::DOFsManager::DOFsData& trial_dofs_data,
                                                                                           const LocalSpace_PCC_2D::ReferenceElement_Data& trial_reference_element_data,
                                                                                           const Eigen::VectorXd& numerical_solution,
                                                                                           const Eigen::VectorXd& numerical_solution_strong,
                                                                                           const std::function<double (const double&, const double&, const double&)> exact_solution_function,
                                                                                           const std::function<std::array<double, 3>(const double &, const double &, const double &)> exact_gradient_solution_function)
        {
          Evaluate_Solution_On_Quadrature_Points_Data result;

          const auto num_solution = to_Eigen_Array(numerical_solution);
          const auto num_solution_strong = to_Eigen_Array(numerical_solution_strong);

          unsigned int num_total_quadrature_points;
          std::vector<unsigned int> cell2D_num_quadrature_points(mesh.Cell2DTotalNumber());
          std::vector<Eigen::MatrixXd> quadrature_points(mesh.Cell2DTotalNumber());
          std::vector<Eigen::VectorXd> quadrature_weigths(mesh.Cell2DTotalNumber());
          std::vector<Eigen::VectorXd> numeric_solution(mesh.Cell2DTotalNumber());
          std::vector<std::array<Eigen::VectorXd, 3>> numeric_gradient_solution(mesh.Cell2DTotalNumber());
          std::vector<Eigen::VectorXd> exact_solution(mesh.Cell2DTotalNumber());
          std::vector<std::array<Eigen::VectorXd, 3>> exact_gradient_solution(mesh.Cell2DTotalNumber());


          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto trial_local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 trial_reference_element_data);

            const auto trial_basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(trial_reference_element_data,
                                                                           trial_local_space_data,
                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0k);

            const auto trial_basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(trial_reference_element_data, trial_local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(trial_reference_element_data, trial_local_space_data);

            const unsigned int num_quadrature_points = cell2D_internal_quadrature.Points.cols();
            num_total_quadrature_points += num_quadrature_points;
            cell2D_num_quadrature_points.at(c) = num_quadrature_points;

            quadrature_points.at(c) = cell2D_internal_quadrature.Points;
            quadrature_weigths.at(c) = cell2D_internal_quadrature.Weights;

            exact_solution.at(c) = function_evaluation(cell2D_internal_quadrature.Points, exact_solution_function);
            exact_gradient_solution.at(c) = function_evaluation(cell2D_internal_quadrature.Points, exact_gradient_solution_function);

            const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, trial_dofs_data);
            const Eigen::VectorXd dofs_values =
                PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                    trial_dofs_data,
                                                                                    local_count_dofs.num_total_dofs,
                                                                                    local_count_dofs.offsets_DOFs,
                                                                                    {0},
                                                                                    {0},
                                                                                    num_solution,
                                                                                    num_solution_strong);

            numeric_solution.at(c) = trial_basis_functions_values * dofs_values;
            numeric_gradient_solution.at(c).at(0) = trial_basis_functions_derivative_values.at(0) * dofs_values;
            numeric_gradient_solution.at(c).at(1) = trial_basis_functions_derivative_values.at(1) * dofs_values;
            numeric_gradient_solution.at(c).at(2) = Eigen::VectorXd::Zero(num_quadrature_points);
          }

          result.quadrature_points.resize(3, num_total_quadrature_points);
          result.quadrature_weigths.resize(num_total_quadrature_points);
          result.numeric_solution.resize(num_total_quadrature_points);
          result.exact_solution.resize(num_total_quadrature_points);
          result.numeric_gradient_solution.at(0).resize(num_total_quadrature_points);
          result.numeric_gradient_solution.at(1).resize(num_total_quadrature_points);
          result.numeric_gradient_solution.at(2).resize(num_total_quadrature_points);
          result.exact_gradient_solution.at(0).resize(num_total_quadrature_points);
          result.exact_gradient_solution.at(1).resize(num_total_quadrature_points);
          result.exact_gradient_solution.at(2).resize(num_total_quadrature_points);

          unsigned int local_num_q = 0;
          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto& cell2D_num_q = cell2D_num_quadrature_points.at(c);

            result.quadrature_points.block(0,
                                           local_num_q,
                                           3,
                                           cell2D_num_q) = quadrature_points.at(c);
            result.quadrature_weigths.segment(local_num_q,
                                              cell2D_num_q) = quadrature_weigths.at(c);
            result.numeric_solution.segment(local_num_q,
                                              cell2D_num_q) = numeric_solution.at(c);
            result.exact_solution.segment(local_num_q,
                                              cell2D_num_q) = exact_solution.at(c);
            result.numeric_gradient_solution.at(0).segment(local_num_q,
                                              cell2D_num_q) = numeric_gradient_solution.at(c).at(0);
            result.numeric_gradient_solution.at(1).segment(local_num_q,
                                              cell2D_num_q) = numeric_gradient_solution.at(c).at(1);
            result.numeric_gradient_solution.at(2).segment(local_num_q,
                                              cell2D_num_q) = numeric_gradient_solution.at(c).at(2);
            result.exact_gradient_solution.at(0).segment(local_num_q,
                                              cell2D_num_q) = exact_gradient_solution.at(c).at(0);
            result.exact_gradient_solution.at(1).segment(local_num_q,
                                              cell2D_num_q) = exact_gradient_solution.at(c).at(1);
            result.exact_gradient_solution.at(2).segment(local_num_q,
                                              cell2D_num_q) = exact_gradient_solution.at(c).at(2);

            local_num_q += cell2D_num_q;
          }

          return result;
        }
        // ***************************************************************************
      } // namespace PCC_2D
    } // namespace Assembler_Utilities
  } // namespace PDETools
} // namespace Polydim
