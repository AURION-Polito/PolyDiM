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

#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
  namespace FEM
  {
    namespace PCC
    {
      // ***************************************************************************
      FEM_Tetrahedron_PCC_3D_LocalSpace_Data FEM_Tetrahedron_PCC_3D_LocalSpace::CreateLocalSpace(
          const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
          const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
      {
        FEM_Tetrahedron_PCC_3D_LocalSpace_Data localSpace;

        Gedim::GeometryUtilitiesConfig geometry_utilities_config;
        geometry_utilities_config.Tolerance1D = polyhedron.Tolerance1D;
        geometry_utilities_config.Tolerance2D = polyhedron.Tolerance2D;
        geometry_utilities_config.Tolerance3D = polyhedron.Tolerance3D;
        Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

        Gedim::MapTetrahedron mapTetrehedron(geometry_utilities);
        localSpace.MapData = mapTetrehedron.Compute(polyhedron.Vertices);

        localSpace.Order = reference_element_data.Order;
        localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

        for (unsigned int e = 0; e < 6; ++e)
        {
          const unsigned int edge_origin_index = polyhedron.Edges(0, e);
          const unsigned int edge_end_index = polyhedron.Edges(1, e);

          const auto &reference_edge = reference_element_data.Edges_by_vertices.at({edge_origin_index, edge_end_index});

          localSpace.polyhedron_to_reference_edge_index[e] = reference_edge.first;
          localSpace.polyhedron_to_reference_edge_direction[e] = reference_edge.second;
        }

        FEM_Triangle_PCC_2D_LocalSpace face_local_space;

        for (unsigned int f = 0; f < 4; ++f)
        {
          const auto &face_geometry = polyhedron.Faces_2D_Geometry[f];

          FEM_PCC_2D_Polygon_Geometry fem_face_geometry = {polyhedron.Tolerance1D,
                                                           polyhedron.Tolerance2D,
                                                           face_geometry.Vertices,
                                                           face_geometry.EdgesDirection,
                                                           face_geometry.EdgesTangent,
                                                           face_geometry.EdgesLength};

          localSpace.Boundary_LocalSpace_Data[f] =
              face_local_space.CreateLocalSpace(reference_element_data.BoundaryReferenceElement_Data, fem_face_geometry);
        }

        for (unsigned int f = 0; f < 4; ++f)
        {
          const unsigned int face_edge_index = polyhedron.Faces[f](1, 0);
          const unsigned int face_edge_next_index = polyhedron.Faces[f](1, 1);
          const unsigned int face_reference_edge_index = localSpace.polyhedron_to_reference_edge_index[face_edge_index];
          const unsigned int face_reference_edge_next_index = localSpace.polyhedron_to_reference_edge_index[face_edge_next_index];
          const unsigned int face_vertex_index = polyhedron.Faces[f](0, 2);

          const auto &ref_face_by_ev =
              reference_element_data.Faces_by_edge_vertex.at({face_reference_edge_index, face_vertex_index});
          const auto &ref_face_by_ee =
              reference_element_data.Faces_by_edges.at({face_reference_edge_index, face_reference_edge_next_index});
          assert(ref_face_by_ee.first.at(0) == ref_face_by_ev);

          localSpace.polyhedron_to_reference_face_index[f] = ref_face_by_ev;
          localSpace.polyhedron_to_reference_face_direction[f] = ref_face_by_ee.second;
          localSpace.polyhedron_to_reference_face_starting_index[f] = ref_face_by_ee.first.at(1);
        }

        localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);

        localSpace.Dof0DsIndex.fill(0);
        for (unsigned int v = 0; v < 4; v++)
        {
          localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

          unsigned int vertex_dofCounter = reference_element_data.NumDofs0D * v;
          for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
          {
            localSpace.DofsMeshOrder[vertex_dofCounter] = d;
            vertex_dofCounter++;
          }
        }

        localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[4]);
        for (unsigned int e = 0; e < 6; ++e)
        {
          localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

          const unsigned int ref_e = localSpace.polyhedron_to_reference_edge_index[e];
          const bool ref_e_direction = localSpace.polyhedron_to_reference_edge_direction[e];
          unsigned int edge_dof_counter = reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * ref_e;

          if (polyhedron.EdgesDirection.at(e) == ref_e_direction)
          {
            for (unsigned int d = localSpace.Dof1DsIndex[e]; d < localSpace.Dof1DsIndex[e + 1]; d++)
            {
              localSpace.DofsMeshOrder[edge_dof_counter] = d;
              edge_dof_counter++;
            }
          }
          else
          {
            for (unsigned int d = localSpace.Dof1DsIndex[e + 1] - 1; d < UINT_MAX && d >= localSpace.Dof1DsIndex[e]; d--)
            {
              localSpace.DofsMeshOrder[edge_dof_counter] = d;
              edge_dof_counter++;
            }
          }
        }

        localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[6]);
        for (unsigned int f = 0; f < 4; ++f)
        {
          localSpace.Dof2DsIndex[f + 1] = localSpace.Dof2DsIndex[f] +
                                          reference_element_data.NumDofs2D;

          const unsigned int ref_f = localSpace.polyhedron_to_reference_face_index[f];
          const bool ref_face_dir = localSpace.polyhedron_to_reference_face_direction[f];
          unsigned int ref_face_s_i = localSpace.polyhedron_to_reference_face_starting_index[f];

          unsigned int face_dof_counter = reference_element_data.NumDofs0D * 4 +
                                          reference_element_data.NumDofs1D * 6 +
                                          reference_element_data.NumDofs2D * ref_f;

          const unsigned int local_space_2D_face_dof_index = 3 * reference_element_data.NumDofs0D +
                                                             3 * reference_element_data.NumDofs1D;

          const auto& local_space_2D = localSpace.Boundary_LocalSpace_Data[f];
          const Eigen::MatrixXd face_internal_dofs_2D = local_space_2D.Dofs.block(0,
                                                                                  local_space_2D_face_dof_index,
                                                                                  3,
                                                                                  reference_element_data.NumDofs2D);
          const auto& face_rotation_matrix = polyhedron.FacesRotationMatrix.at(f);
          const auto& face_translation = polyhedron.FacesTranslation.at(f);
          const auto face_internal_dofs_3D_from_reference_2D = geometry_utilities.RotatePointsFrom2DTo3D(face_internal_dofs_2D,
                                                                                                         face_rotation_matrix,
                                                                                                         face_translation);



          const Eigen::MatrixXd reference_face_dofs = reference_element_data.DofPositions.block(0,
                                                                                                face_dof_counter,
                                                                                                3,
                                                                                                reference_element_data.NumDofs2D);


          const auto face_internal_dofs_3D_from_reference_3D = Gedim::MapTetrahedron::F(localSpace.MapData, reference_face_dofs);


          std::cout << "\tFace " << f << " dir " << polyhedron.FacesDirection.at(f) << " dir " << ref_face_dir << " s_i " << ref_face_s_i << " ";
          std::cout << "Ref F " << ref_f<< " ";
          std::cout << "vertices " << polyhedron.Faces[f](0, 0) << ", ";
          std::cout << polyhedron.Faces[f](0, 1) << ", ";
          std::cout << polyhedron.Faces[f](0, 2) << " ";
          std::cout << "edges " << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 0)] << ", ";
          std::cout << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 1)] << ", ";
          std::cout << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 2)]<< std::endl;

//          std::cout<< "inter_from_2D:\n"<< face_internal_dofs_3D_from_reference_2D<< std::endl;
//          std::cout<< "inter_from_3D:\n"<< face_internal_dofs_3D_from_reference_3D<< std::endl;

          std::array<unsigned int, 3> original_dmo;
          unsigned int s_h = 0;
          for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
          {
            original_dmo[face_dof_counter + s_h] = d;
            s_h++;
          }
          std::cout << " o_dmo " << original_dmo[face_dof_counter] << ", ";
          std::cout << original_dmo[face_dof_counter + 1] << ", ";
          std::cout << original_dmo[face_dof_counter + 2] << " ";

          unsigned int shift = 0;
          for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
          {
            const auto find_point_index = geometry_utilities.FindPointInPoints(face_internal_dofs_3D_from_reference_3D,
                                                                               face_internal_dofs_3D_from_reference_2D.col(shift));
            assert(find_point_index.size() == 1);

            const unsigned int find_shift = find_point_index[0];
            localSpace.DofsMeshOrder[face_dof_counter + find_shift] = d;
            shift++;
          }

          std::cout << " dmo " << localSpace.DofsMeshOrder[face_dof_counter] << ", ";
          std::cout << localSpace.DofsMeshOrder[face_dof_counter + 1] << ", ";
          std::cout << localSpace.DofsMeshOrder[face_dof_counter + 2] << std::endl;

        }

        /* TOREMOVE
    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[6]);
    for (unsigned int f = 0; f < 4; ++f)
    {
        localSpace.Dof2DsIndex[f + 1] = localSpace.Dof2DsIndex[f] + reference_element_data.NumDofs2D;

        const unsigned int ref_f = localSpace.polyhedron_to_reference_face_index[f];
        const bool ref_face_dir = localSpace.polyhedron_to_reference_face_direction[f];
        unsigned int ref_face_s_i = localSpace.polyhedron_to_reference_face_starting_index[f];
        unsigned int face_dof_counter = reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * 6 +
                                        reference_element_data.NumDofs2D * ref_f;

        // std::cout << "\tFace " << f << " dir " << polyhedron.FacesDirection.at(f) << " ";
        // std::cout << "Ref F " << ref_f << " dir " << ref_face_dir << " s_i " << ref_face_s_i << " ";
        // std::cout << "vertices " << polyhedron.Faces[f](0, 0) << ", ";
        // std::cout << polyhedron.Faces[f](0, 1) << ", ";
        // std::cout << polyhedron.Faces[f](0, 2) << " ";
        // std::cout << "edges " << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 0)] << ", ";
        // std::cout << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 1)] << ", ";
        // std::cout << localSpace.polyhedron_to_reference_edge_index[polyhedron.Faces[f](1, 2)];

        const std::array<unsigned int, 3> shift_map = {0, 2, 1};
        const std::array<unsigned int, 3> shift_map_rev = {1, 2, 0};
        unsigned int shift = ref_face_dir ? shift_map[ref_face_s_i] : shift_map_rev[ref_face_s_i];

        if (ref_face_dir)
        {
            if (ref_f == 3)
            {
                for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
                {
                    unsigned int loc_shift = (shift + 2) % 3;
                    localSpace.DofsMeshOrder[face_dof_counter + loc_shift] = d;
                    shift++;
                }
            }
            else
            {
                for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
                {
                    unsigned int loc_shift = shift % 3;
                    localSpace.DofsMeshOrder[face_dof_counter + loc_shift] = d;
                    shift++;
                }
            }
        }
        else
        {
            if (ref_f == 3)
            {
                for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f]; d--)
                {
                    unsigned int loc_shift = (shift + 2) % 3;
                    localSpace.DofsMeshOrder[face_dof_counter + loc_shift] = d;
                    shift++;
                }
            }
            else
            {
                for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f]; d--)
                {
                    unsigned int loc_shift = shift % 3;
                    localSpace.DofsMeshOrder[face_dof_counter + loc_shift] = d;
                    shift++;
                }
            }
        }

        std::array<unsigned int, 3> original_dmo;
        unsigned int s_h = 0;
        for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
        {
            original_dmo[face_dof_counter + s_h] = d;
            s_h++;
        }

        // std::cout << " o_dmo " << original_dmo[face_dof_counter] << ", ";
        // std::cout << original_dmo[face_dof_counter + 1] << ", ";
        // std::cout << original_dmo[face_dof_counter + 2] << " ";

        // std::cout << " dmo " << localSpace.DofsMeshOrder[face_dof_counter] << ", ";
        // std::cout << localSpace.DofsMeshOrder[face_dof_counter + 1] << ", ";
        // std::cout << localSpace.DofsMeshOrder[face_dof_counter + 2] << std::endl;

        // if (polyhedron.FacesDirection.at(f))
        //{
        //   for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
        //   {
        //       localSpace.DofsMeshOrder[face_dof_counter] = d;
        //       face_dof_counter++;
        //   }
        // }
        // else
        //{
        //     for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f];
        //     d--)
        //     {
        //         localSpace.DofsMeshOrder[face_dof_counter] = d;
        //         face_dof_counter++;
        //     }
        // }
    }
    TO REMOVE */

        localSpace.Dof3DsIndex.fill(localSpace.Dof2DsIndex[4]);
        localSpace.Dof3DsIndex[1] = localSpace.Dof3DsIndex[0] + reference_element_data.NumDofs3D;

        unsigned int cell_dof_counter =
            reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * 6 + reference_element_data.NumDofs2D * 4;
        for (unsigned int d = localSpace.Dof3DsIndex[0]; d < localSpace.Dof3DsIndex[1]; d++)
        {
          localSpace.DofsMeshOrder[cell_dof_counter] = d;
          cell_dof_counter++;
        }

        localSpace.Dofs = MapValues(localSpace, Gedim::MapTetrahedron::F(localSpace.MapData, reference_element_data.DofPositions));

        localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTetrahedronQuadrature, localSpace.MapData);
        localSpace.BoundaryQuadrature = BoundaryQuadrature(localSpace.Boundary_LocalSpace_Data, polyhedron);

        return localSpace;
      }
      // ***************************************************************************
      MatrixXd FEM_Tetrahedron_PCC_3D_LocalSpace::MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                            const Eigen::MatrixXd &referenceValues) const
      {
        Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

        for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
          basisFunctionValuesOrdered.col(local_space.DofsMeshOrder.at(d)) << referenceValues.col(d);

        return basisFunctionValuesOrdered;
      }
      // ***************************************************************************
      std::vector<MatrixXd> FEM_Tetrahedron_PCC_3D_LocalSpace::MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                                   const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
      {
        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
              3,
              Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

        for (unsigned int i = 0; i < 3; i++)
        {
          basisFunctionsDerivativeValues[i] = local_space.MapData.QInv(i, i) * referenceDerivateValues[i];
          for (unsigned int j = 0; j < i; j++)
          {
            basisFunctionsDerivativeValues[i] += local_space.MapData.QInv(j, i) * referenceDerivateValues[j];
            basisFunctionsDerivativeValues[j] += local_space.MapData.QInv(i, j) * referenceDerivateValues[i];
          }
        }

        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
              3,
              Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

        for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
        {
          const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
          basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
          basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
          basisFunctionsDerivativeValuesOrdered.at(2).col(dofOrder) << basisFunctionsDerivativeValues.at(2).col(d);
        }

        return basisFunctionsDerivativeValuesOrdered;
      }
      // ***************************************************************************
      Gedim::Quadrature::QuadratureData FEM_Tetrahedron_PCC_3D_LocalSpace::InternalQuadrature(
          const Gedim::Quadrature::QuadratureData &reference_quadrature,
          const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const
      {
        Gedim::Quadrature::QuadratureData quadrature;

        quadrature.Points = Gedim::MapTetrahedron::F(mapData, reference_quadrature.Points);
        quadrature.Weights = reference_quadrature.Weights.array() *
                             Gedim::MapTetrahedron::DetJ(mapData, reference_quadrature.Points).array().abs();

        return quadrature;
      }
      // ***************************************************************************
      std::array<Gedim::Quadrature::QuadratureData, 4> FEM_Tetrahedron_PCC_3D_LocalSpace::BoundaryQuadrature(
          const std::array<FEM_Triangle_PCC_2D_LocalSpace_Data, 4> &faces_local_space_data,
          const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
      {
        std::array<Gedim::Quadrature::QuadratureData, 4> faces_quadrature;

        for (unsigned int f = 0; f < 4; ++f)
        {
          auto &face_quadrature = faces_quadrature.at(f);
          face_quadrature = faces_local_space_data[f].InternalQuadrature;

          const Eigen::Vector3d &face_translation = polyhedron.FacesTranslation[f];
          const Eigen::Matrix3d &face_rotation = polyhedron.FacesRotationMatrix[f];

          face_quadrature.Points = (face_rotation * face_quadrature.Points).colwise() + face_translation;
        }

        return faces_quadrature;
      }
      // ***************************************************************************
    } // namespace PCC

  } // namespace FEM
} // namespace Polydim
