#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"

namespace Elliptic_MCC_2D
{
struct Program_configuration final
{
    enum struct MeshGenerators
    {
        Tri = 0, // triangular mesh
        OFFImporter = 1, // imported off mesh
    };

    enum struct VemTypes
    {
        Vem = 1,
        VemPartial = 2,
        VemOrtho = 3
    };

    enum struct ProgramTypes
    {
        Poisson = 0,
        PatchTest = 1
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
                                           "Mesh 2D gereator type, 0 - triangle; 1 - OFF Importer; (Default: 0)");

        Gedim::Configurations::AddProperty("MeshOFF_FilePath",
                                           "./",
                                           "Mesh OFF imported file path, use it when meshIsConcave=true (Default: './')");

        Gedim::Configurations::AddProperty("MeshMaxArea",
                                           0.1,
                                           "Mesh 2D maximum cell area (Default: 0.1)");

        /// Method parameters
        Gedim::Configurations::AddProperty("VemOrder",
                                           static_cast<unsigned int>(0),
                                           "VEM order (Default: 0)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance",
                                           true,
                                           "Compute VEM Performance (Default: true)");
        Gedim::Configurations::AddProperty("ComputeConditionNumber",
                                           false,
                                           "Compute Condition Number (Default: false)");


        /// Program parameters
        Gedim::Configurations::AddProperty("VemType",
                                           static_cast<unsigned int>(VemTypes::Vem),
                                           "VEM type, 1 - Vem; 2 - VemPartial; 3 - VemOrtho (Default: 1)");
        Gedim::Configurations::AddProperty("ProgramType",
                                           static_cast<unsigned int>(ProgramTypes::Poisson),
                                           "Program type, 0 - Poisson; 1 - PatchTest; (Default: 0)");
    }

    inline string ExportFolder() const
    { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }

    inline double GeometricTolerance() const
    { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }

    inline MeshGenerators MeshGenerator() const
    { return (MeshGenerators)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"); }
    inline std::string MeshOFF_FilePath() const
    { return Gedim::Configurations::GetPropertyValue<string>("MeshOFF_FilePath"); }
    inline double MeshMaxArea() const
    { return Gedim::Configurations::GetPropertyValue<double>("MeshMaxArea"); }

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
