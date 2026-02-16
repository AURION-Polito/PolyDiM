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
                                              const DOFs::DOFsManager::DOFsData &dofs_data,
                                              const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                              const std::function<double(const double &, const double &, const double &, const Eigen::VectorXd &)> source_term_function)
        {
          Gedim::Eigen_Array<> forcing_term;

          forcing_term.SetSize(dofs_data.NumberDOFs);

          Polydim::PDETools::Equations::EllipticEquation equation;

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            const Eigen::VectorXd source_term_values = function_evaluation(cell2D_internal_quadrature.Points, source_term_function);

            const Eigen::VectorXd local_rhs =
                equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

            assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) ==
                   dofs_data.CellsGlobalDOFs[2].at(c).size());

            Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                              local_matrix_to_global_matrix_dofs_data,
                                                                                              local_rhs,
                                                                                              forcing_term);
          }

          forcing_term.Create();

          return static_cast<Eigen::VectorXd &>(forcing_term);
        }
        // ***************************************************************************
        Variational_Operator assembler_elliptic_operator(
            const Gedim::GeometryUtilities &geometry_utilities,
            const Gedim::MeshMatricesDAO &mesh,
            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
            const DOFs::DOFsManager::DOFsData &dofs_data,
            const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
            const std::function<double(const double &, const double &, const double &, const Eigen::VectorXd &)> diffusion_term_function)
        {
          Gedim::Eigen_SparseArray<> elliptic_matrix;
          Gedim::Eigen_SparseArray<> elliptic_strong_matrix;

          elliptic_matrix.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::None);
          elliptic_strong_matrix.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);

          Polydim::PDETools::Equations::EllipticEquation equation;

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

            const auto basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            const Eigen::VectorXd diffusion_term_values = function_evaluation(cell2D_internal_quadrature.Points, diffusion_term_function);

            const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                                basis_functions_derivative_values,
                                                                                cell2D_internal_quadrature.Weights);

            const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
            const Eigen::MatrixXd local_A_stab =
                k_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data, local_space_data);

            assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) ==
                   dofs_data.CellsGlobalDOFs[2].at(c).size());

            Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                              local_matrix_to_global_matrix_dofs_data,
                                                                                              local_matrix_to_global_matrix_dofs_data,
                                                                                              local_A + local_A_stab,
                                                                                              elliptic_matrix,
                                                                                              elliptic_strong_matrix);
          }

          elliptic_matrix.Create();
          elliptic_strong_matrix.Create();

          return convert_operator(elliptic_matrix, elliptic_strong_matrix);
        }
        // ***************************************************************************
        Eigen::VectorXd assembler_strong_solution(const Gedim::GeometryUtilities &geometry_utilities,
                                                  const Gedim::MeshMatricesDAO &mesh,
                                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                  const DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                                  const DOFs::DOFsManager::DOFsData &dofs_data,
                                                  const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                  const std::function<double(const unsigned int, const double &, const double &, const double &)> strong_solution_function)
        {
          Gedim::Eigen_Array<> strong_solution;

          strong_solution.SetSize(dofs_data.NumberStrongs);

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) ==
                   dofs_data.CellsGlobalDOFs[2].at(c).size());

            for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); ++v)
            {
              const unsigned int cell0D_index = mesh.Cell2DVertex(c, v);
              const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                continue;

              const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

              const auto strong_boundary_values = function_evaluation(boundary_info.Marker, coordinates, strong_solution_function);

              const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

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

            for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
            {
              const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

              const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
              const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
                  local_dofs.size() == 0)
                continue;

              const auto edge_dofs_coordinates =
                  Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

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

          strong_solution.Create();

          return static_cast<Eigen::VectorXd &>(strong_solution);
        }
        // ***************************************************************************
        Exact_Solution_Data assembler_exact_solution(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const std::function<double(const double &, const double &, const double &)> exact_solution_function)
        {
          Gedim::Eigen_Array<> exact_solution;
          Gedim::Eigen_Array<> exact_solution_strong;

          exact_solution.SetSize(dofs_data.NumberDOFs);
          exact_solution_strong.SetSize(dofs_data.NumberStrongs);

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

              const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);
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
            if (reference_element_data.Order > 1)
            {
              const auto local_space_data =
                  Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                         geometry_utilities.Tolerance2D(),
                                                                         mesh_geometric_data,
                                                                         c,
                                                                         reference_element_data);

              // Assemble strong boundary condition on Cell1Ds
              for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
              {
                const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

                const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

                const auto edge_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

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

              const auto local_dofs = dofs_data.CellsDOFs.at(2).at(c);

              if (local_dofs.size())
              {
                const auto internal_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofsCoordinates(reference_element_data, local_space_data);

                const Eigen::VectorXd initial_values_at_dofs =
                    function_evaluation(internal_dofs_coordinates.Points, exact_solution_function);

                const Eigen::VectorXd dofs_internal =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofs(reference_element_data,
                                                                       local_space_data,
                                                                       initial_values_at_dofs,
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
        Eigen::VectorXd assembler_weak_term(const Gedim::GeometryUtilities &geometry_utilities,
                                            const Gedim::MeshMatricesDAO &mesh,
                                            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                            const DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                            const DOFs::DOFsManager::DOFsData &dofs_data,
                                            const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                            const std::function<double(const unsigned int, const double &, const double &, const double &)> weak_term_function)
        {
          Gedim::Eigen_Array<> weak_term;

          weak_term.SetSize(dofs_data.NumberDOFs);

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
          {{std::cref(dofs_data)}, {0}, {0}, {0}};

          for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
          {
            if (!mesh.Cell2DIsActive(c))
              continue;

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(c).cols();

            for (unsigned int ed = 0; ed < numVertices; ed++)
            {
              const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

              const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
                continue;

              // compute vem values
              const auto weakReferenceSegment =
                  Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

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
                  Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValuesOnEdge(ed, reference_element_data, local_space_data, pointsCurvilinearCoordinates);

              // compute values of Neumann condition
              const Eigen::VectorXd neumannContributions =
                  weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues;

              for (unsigned int p = 0; p < 2; ++p)
              {
                const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index, p);

                const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

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

              const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
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
        Post_Process_Data_Cell0Ds assembler_extract_cell0Ds(const Gedim::MeshMatricesDAO &mesh,
                                                            const DOFs::DOFsManager::DOFsData &dofs_data,
                                                            const Eigen::VectorXd &numerical_solution,
                                                            const Eigen::VectorXd &numerical_solution_strong,
                                                            const std::function<double(const double &, const double &, const double &)> exact_solution_function)
        {
          Post_Process_Data_Cell0Ds result;

          const auto num_solution = to_Eigen_Array(numerical_solution);
          const auto num_solution_strong = to_Eigen_Array(numerical_solution_strong);

          result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
          result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

          for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
          {
            if (!mesh.Cell0DIsActive(p))
              continue;

            result.cell0Ds_exact[p] = function_evaluation(mesh.Cell0DCoordinates(p), exact_solution_function)[0];

            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
              const auto &local_dof_i = local_dofs.at(loc_i);

              switch (local_dof_i.Type)
              {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                  result.cell0Ds_numeric[p] = numerical_solution_strong[local_dof_i.Global_Index];
                  break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                  result.cell0Ds_numeric[p] = numerical_solution[local_dof_i.Global_Index];
                  break;
                default:
                  throw std::runtime_error("Unknown DOF Type");
              }
            }
          }

          return result;
        }
        // ***************************************************************************
        Post_Process_Data_ErrorL2 assembler_error_L2(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
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

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data,
                                                                           local_space_data,
                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0k);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            const auto exact_solution_values = function_evaluation(cell2D_internal_quadrature.Points, exact_solution_function);

            const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
            const Eigen::VectorXd dofs_values =
                PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                    dofs_data,
                                                                                    local_count_dofs.num_total_dofs,
                                                                                    local_count_dofs.offsets_DOFs,
                                                                                    {0},
                                                                                    {0},
                                                                                    num_solution,
                                                                                    num_solution_strong);

            const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
            const Eigen::VectorXd local_numeric_norm_L2 = (basis_functions_values * dofs_values).array().square();
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
        Post_Process_Data_ErrorH1 assembler_error_H1(const Gedim::GeometryUtilities &geometry_utilities,
                                                     const Gedim::MeshMatricesDAO &mesh,
                                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                     const DOFs::DOFsManager::DOFsData &dofs_data,
                                                     const LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
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

            const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                                 geometry_utilities.Tolerance2D(),
                                                                                                 mesh_geometric_data,
                                                                                                 c,
                                                                                                 reference_element_data);

            const auto basis_functions_derivative_values =
                Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

            const auto cell2D_internal_quadrature =
                Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

            const auto exact_derivative_solution_values = function_evaluation(cell2D_internal_quadrature.Points, exact_gradient_solution_function);

            const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
            const Eigen::VectorXd dofs_values =
                PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                    dofs_data,
                                                                                    local_count_dofs.num_total_dofs,
                                                                                    local_count_dofs.offsets_DOFs,
                                                                                    {0},
                                                                                    {0},
                                                                                    num_solution,
                                                                                    num_solution_strong);

            const Eigen::VectorXd local_error_H1 =
                (basis_functions_derivative_values.at(0) * dofs_values -
                exact_derivative_solution_values.at(0)).array().square() +
                (basis_functions_derivative_values.at(1) * dofs_values -
                exact_derivative_solution_values.at(1)).array().square();

            const Eigen::VectorXd local_numeric_norm_H1 =
                (basis_functions_derivative_values.at(0) * dofs_values).array().square() +
                (basis_functions_derivative_values.at(1) * dofs_values).array().square();

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
      } // namespace PCC_2D
    } // namespace Assembler_Utilities
  } // namespace PDETools
} // namespace Polydim
