#ifndef __PDETOOLS_MESH_PDE_Mesh_Generation_HPP
#define __PDETOOLS_MESH_PDE_Mesh_Generation_HPP

#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

namespace Polydim
{
  namespace PDETools
  {
    namespace Mesh
    {
      namespace PDE_Mesh_Generation
      {
        struct PDE_Domain_2D final
        {
            Eigen::MatrixXd vertices;
            double area;
        };

        enum struct MeshGenerator_Types_2D
        {
          Triangular = 0, ///< generated triangular mesh
          Minimal = 1, ///< generated minimal mesh
          Polygonal = 2, ///< generated voronoi polygonal mesh
          OFFImporter = 3 ///< imported off mesh
        };

        void create_mesh_2D(const Gedim::GeometryUtilities& geometry_utilities,
                            const Gedim::MeshUtilities& mesh_utilities,
                            const MeshGenerator_Types_2D& mesh_type,
                            const PDE_Domain_2D& pde_domain,
                            const double& max_relative_area,
                            Gedim::MeshMatricesDAO& mesh)
        {
          switch (mesh_type)
          {
            case MeshGenerator_Types_2D::Triangular:
            {
              const double max_cell_area = pde_domain.area *
                                           max_relative_area;
              mesh_utilities.CreateTriangularMesh(pde_domain.vertices,
                                                  max_cell_area,
                                                  mesh);
            }
              break;
            case MeshGenerator_Types_2D::Minimal:
            {
              mesh_utilities.Mesh2DFromPolygon(pde_domain.vertices,
                                               {},
                                               {},
                                               mesh);
            }
              break;
            case MeshGenerator_Types_2D::Polygonal:
            {
              const unsigned num_cells =
                  static_cast<unsigned int>(std::max(1.0, 1.0 /
                                                     max_relative_area));


              mesh_utilities.CreatePolygonalMesh(geometry_utilities,
                                                 pde_domain.vertices,
                                                 num_cells,
                                                 10,
                                                 mesh);
            }
              break;
            default:
              throw std::runtime_error("MeshGenerator_Types_2D " +
                                       std::to_string((unsigned int)mesh_type) +
                                       " not supported");
          }
        }

        void import_mesh_2D(const Gedim::GeometryUtilities& geometry_utilities,
                            const Gedim::MeshUtilities& mesh_utilities,
                            const MeshGenerator_Types_2D& mesh_type,
                            const std::string& file_path,
                            Gedim::MeshMatricesDAO& mesh)
        {
          switch (mesh_type)
          {
            case MeshGenerator_Types_2D::OFFImporter:
            {
              mesh_utilities.ImportObjectFileFormat(file_path,
                                                    mesh);
            }
              break;
            default:
              throw std::runtime_error("MeshGenerator_Types_2D " +
                                       std::to_string((unsigned int)mesh_type) +
                                       " not supported");
          }
        }
      }
    }
  }
}

#endif
