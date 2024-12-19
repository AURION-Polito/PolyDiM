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
    Gedim::MeshMatricesDAO &mesh_data;

    inline unsigned int Cell0Ds_number() const
    {
        return mesh_data.Cell0DTotalNumber();
    }
    inline unsigned int Cell1Ds_number() const
    {
        return mesh_data.Cell1DTotalNumber();
    }
    inline unsigned int Cell2Ds_number() const
    {
        return mesh_data.Cell2DTotalNumber();
    }
    inline unsigned int Cell3Ds_number() const
    {
        return mesh_data.Cell3DTotalNumber();
    }

    inline unsigned int Cell0D_marker(const unsigned int cell0D_index) const
    {
        return mesh_data.Cell0DMarker(cell0D_index);
    }
    inline unsigned int Cell1D_marker(const unsigned int cell1D_index) const
    {
        return mesh_data.Cell1DMarker(cell1D_index);
    }
    inline unsigned int Cell2D_marker(const unsigned int cell2D_index) const
    {
        return mesh_data.Cell2DMarker(cell2D_index);
    }
    inline unsigned int Cell3D_marker(const unsigned int cell3D_index) const
    {
        return mesh_data.Cell3DMarker(cell3D_index);
    }

    inline std::array<unsigned int, 2> Cell1D_vertices(const unsigned int cell1D_index) const
    {
        return {mesh_data.Cell1DOrigin(cell1D_index), mesh_data.Cell1DEnd(cell1D_index)};
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
        return mesh_data.Cell3DFaces(cell3D_index);
    }
};
} // namespace Mesh
} // namespace PDETools
} // namespace Polydim

#endif
