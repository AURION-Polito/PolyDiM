#ifndef __PDETOOLS_MESH_MeshMatricesDAO_mesh_connectivity_data_HPP
#define __PDETOOLS_MESH_MeshMatricesDAO_mesh_connectivity_data_HPP

#include "MeshMatricesDAO.hpp"

namespace Polydim
{
  namespace PDETools
  {
    namespace Mesh
    {
      struct MeshMatricesDAO_mesh_connectivity_data final
      {
          Gedim::MeshMatricesDAO& mesh_data;

          inline std::array<unsigned int, 2> Cell1D_vertices(const unsigned int cell1D_index) const
          {
            return
            {
              mesh_data.Cell1DOrigin(cell1D_index), mesh_data.Cell1DEnd(cell1D_index)
            };
          }

          inline std::vector<unsigned int> Cell2D_vertices(const unsigned int cell2D_index) const
          {
            return mesh_data.Cell2DVertices(cell2D_index);
          }

          inline std::vector<unsigned int> Cell2D_edges(const unsigned int cell2D_index) const
          {
            return mesh_data.Cell2DEdges(cell2D_index);
          }

          inline std::vector<unsigned int> Cell3D_vertices(const unsigned int cell3D_index) const
          {
            return mesh_data.Cell3DVertices(cell3D_index);
          }

          inline std::vector<unsigned int> Cell3D_edges(const unsigned int cell3D_index) const
          {
            return mesh_data.Cell3DEdges(cell3D_index);
          }

          inline std::vector<unsigned int> Cell3D_faces(const unsigned int cell3D_index) const
          {
            return mesh_data.Cell3DEdges(cell3D_index);
          }
      };
    }
  }
}

#endif
