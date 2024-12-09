#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"

namespace Elliptic_PCC_2D
{
  struct Program_configuration final
  {
    public:
      enum struct MeshGenerators
      {
        Tri = 0, // triangular mesh
        Rect = 1, // rectangular mesh
        CsvImporter = 2, // imported csv mesh
        RectWithHangs = 3, // rectangular mesh with hanging nodes
        RandomVoronoi = 4,
        OFFImporter = 5, // imported off mesh
      };

      enum struct VemTypes
      {
        EVem = 1,
        E2Vem = 2,
        EVemOrtho = 3,
        E2VemOrtho = 4
      };

      enum struct ProgramTypes
      {
        Poisson = 0,
        SinSinEnlarged = 1
      };

      Program_configuration()
      {
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder",
                                           "./Run",
                                           "Folder where to export data (Default: ./Export)");
        // Geometric parameters
        Gedim::Configurations::AddProperty("GeometricTolerance",
                                           1.0e-8,
                                           "Geometric tolerance to perform 1D operations (Default: machine epsilon)");
        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshGenerator",
                                           static_cast<unsigned int>(MeshGenerators::Tri),
                                           "Mesh 2D gereator type, 0 - triangle; 1 - rectangle; 2 - CSV Importer; 3 - rectangle plus N nodes; 4 - voronoi; 5 - OFF Importer; (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImporterFolder",
                                           "./",
                                           "Mesh importer folder (Default: './')");
        Gedim::Configurations::AddProperty("MeshImporterSeparator",
                                           ';',
                                           "Mesh importer separator (Default: ';')");
        Gedim::Configurations::AddProperty("MeshIsConcave",
                                           false,
                                           "True if the mesh is concave (Default: false)");
        Gedim::Configurations::AddProperty("MeshConvexImporterFolder",
                                           "./",
                                           "Convex Mesh importer folder, use it when meshIsConcave=true (Default: './')");


        Gedim::Configurations::AddProperty("MeshOFF_Original_FilePath",
                                           "./",
                                           "Original Mesh OFF imported file path, use it when meshIsConcave=true (Default: './')");
        Gedim::Configurations::AddProperty("MeshOFF_Aggregated_FilePath",
                                           "./",
                                           "Aggregated Mesh OFF imported file path, use it when meshIsConcave=true (Default: './')");
        Gedim::Configurations::AddProperty("HierarchyMap_FilePath",
                                           "./",
                                           "HierarchyMap OFF imported file path, use it when meshIsConcave=true (Default: './')");



        Gedim::Configurations::AddProperty("MeshMinimumCellSize",
                                           0.1,
                                           "Mesh 2D minimum cell 2D size (Default: 0.1)");

        Gedim::Configurations::AddProperty("MeshNumRectanglesTangentBasis",
                                           static_cast<unsigned int>(1),
                                           "Mesh 2D Number of rectangles tangent to basis (Default: 1)");
        Gedim::Configurations::AddProperty("MeshNumRectanglesTangentHeight",
                                           static_cast<unsigned int>(1),
                                           "Mesh 2D Number of rectangles tangent to height (Default: 1)");

        Gedim::Configurations::AddProperty("MeshNumAddedHangingNodesToEachRectangle",
                                           static_cast<unsigned int>(0),
                                           "Mesh 2D Number of hanging nodes which must be added to each rectangle (Default: 0)");

        Gedim::Configurations::AddProperty("MeshNumPointsVoronoi",
                                           static_cast<unsigned int>(5),
                                           "Mesh 2D Number of points to generate a random Voronoi MEsh 2D (Default: 5)");

        Gedim::Configurations::AddProperty("MeshNumIterationsVoronoi",
                                           static_cast<unsigned int>(2),
                                           "Mesh 2D Number of iterations to generate a random Voronoi MEsh 2D (Default: 2)");




        /// Method parameters
        Gedim::Configurations::AddProperty("VemOrder",
                                           static_cast<unsigned int>(1),
                                           "VEM order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance",
                                           true,
                                           "Compute VEM Performance (Default: true)");
        Gedim::Configurations::AddProperty("ComputeConditionNumber",
                                           false,
                                           "Compute Condition Number (Default: false)");


        /// Program parameters
        Gedim::Configurations::AddProperty("VemType",
                                           static_cast<unsigned int>(VemTypes::EVem),
                                           "VEM type, 1 - EVem; 2 - E2Vem; 3 - EVemOrtho; 4 - E2VemOrtho (Default: 1)");
        Gedim::Configurations::AddProperty("ProgramType",
                                           static_cast<unsigned int>(ProgramTypes::Poisson),
                                           "Program type, 0 - Poisson; 1 - SinSinEnlarged; (Default: 0)");
      }

      inline string ExportFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }

      inline double GeometricTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }

      inline MeshGenerators MeshGenerator() const
      { return (MeshGenerators)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"); }
      inline string MeshImporterFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("MeshImporterFolder"); }
      inline char MeshImporterSeparator() const
      { return Gedim::Configurations::GetPropertyValue<char>("MeshImporterSeparator"); }
      inline bool MeshIsConcave() const
      { return Gedim::Configurations::GetPropertyValue<bool>("MeshIsConcave"); }
      inline string MeshConvexImporterFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("MeshConvexImporterFolder"); }

      inline std::string MeshOFF_Original_FilePath() const
      { return Gedim::Configurations::GetPropertyValue<string>("MeshOFF_Original_FilePath"); }
      inline std::string MeshOFF_Aggregated_FilePath() const
      { return Gedim::Configurations::GetPropertyValue<string>("MeshOFF_Aggregated_FilePath"); }
      inline std::string HierarchyMap_FilePath() const
      { return Gedim::Configurations::GetPropertyValue<string>("HierarchyMap_FilePath"); }

      inline double MeshMinimumCellSize() const
      { return Gedim::Configurations::GetPropertyValue<double>("MeshMinimumCellSize"); }

      inline unsigned int MeshNumRectanglesTangentBasis() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("MeshNumRectanglesTangentBasis"); }

      inline unsigned int MeshNumRectanglesTangentHeight() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("MeshNumRectanglesTangentHeight"); }

      inline unsigned int MeshNumAddedHangingNodesToEachRectangle() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("MeshNumAddedHangingNodesToEachRectangle"); }


      inline unsigned int MeshNumPointsVoronoi() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("MeshNumPointsVoronoi"); }
      inline unsigned int MeshNumIterationsVoronoi() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("MeshNumIterationsVoronoi"); }

      inline bool ComputeConditionNumber() const
      { return Gedim::Configurations::GetPropertyValue<bool>("ComputeConditionNumber"); }
      inline bool ComputeVEMPerformance() const
      { return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance"); }
      inline unsigned int VemOrder() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder"); }

      inline VemTypes VemType() const
      { return (VemTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("VemType"); }

      inline ProgramTypes ProgramType() const
      { return (ProgramTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("ProgramType"); }
  };
}

#endif
