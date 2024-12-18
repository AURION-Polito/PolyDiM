#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
struct Program_configuration final
{
public:
    Program_configuration()
    {
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder",
                                           "./Run",
                                           "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshGenerator",
                                           static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular),
                                           "Mesh 2D gereator type, 0 - Triangular; 1 - Minimal; 2 - Polygonal; 3 - OFF Importer (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath",
                                           "./",
                                           "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxArea",
                                           0.1,
                                           "Mesh 2D maximum relative cell area (Default: 0.1)");

        /// Method parameters
        Gedim::Configurations::AddProperty("VemOrder",
                                           static_cast<unsigned int>(1),
                                           "VEM order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance",
                                           true,
                                           "Compute VEM Performance (Default: true)");

    }

    inline string ExportFolder() const
    { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }

    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D MeshGenerator() const
    { return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"); }
    inline std::string MeshImportFilePath() const
    { return Gedim::Configurations::GetPropertyValue<string>("MeshImportFilePath"); }
    inline double MeshMaxArea() const
    { return Gedim::Configurations::GetPropertyValue<double>("MeshMaxArea"); }

    inline bool ComputeVEMPerformance() const
    { return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance"); }
    inline unsigned int VemOrder() const
    { return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder"); }

};
}
}
}

#endif
