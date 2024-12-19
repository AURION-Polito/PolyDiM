#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "test_definition.hpp"
#include "VEM_PCC_3D_Creator.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
struct Program_configuration final
{
public:
    Program_configuration()
    {
        Gedim::Configurations::AddProperty("TestType",
                                           static_cast<unsigned int>(Polydim::examples::Elliptic_PCC_3D::test::Test_Types::Patch_Test),
                                           "Test Type 1 - Patch_Test; 2 - Poisson_Polynomial_Problem (Default: 1)");
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder",
                                           "./Run",
                                           "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshGenerator",
                                           static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Tetrahedral),
                                           "Mesh 3D gereator type, 0 - Tetrahedral; 1 - Minimal; 2 - Polyhedral; 3 - OVMImporter; 4 - VtkImporter (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath",
                                           "./",
                                           "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxVolume",
                                           0.1,
                                           "Mesh 3D maximum relative cell volume (Default: 0.1)");


        Gedim::Configurations::AddProperty("GeometricTolerance1D",
                                           1.0e-12,
                                           "Geometric Tolerance 1D (Default: 1.0e-12)");

        Gedim::Configurations::AddProperty("GeometricTolerance2D",
                                           1.0e-14,
                                           "Geometric Tolerance 2D (Default: 1.0e-14)");

        Gedim::Configurations::AddProperty("GeometricTolerance3D",
                                           1.0e-15,
                                           "Geometric Tolerance 3D (Default: 1.0e-15)");

        // Method parameters
        Gedim::Configurations::AddProperty("VemType",
                                           static_cast<unsigned int>(Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Types::VEM_PCC_3D_LocalSpace),
                                           "Vem Type, 1 - EVem; 2 - EVem_Inertia; 3 - EVem_Ortho (Default: 1)");
        Gedim::Configurations::AddProperty("VemOrder",
                                           static_cast<unsigned int>(1),
                                           "VEM order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance",
                                           true,
                                           "Compute VEM Performance (Default: true)");

    }

    inline string ExportFolder() const
    { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }

    inline Polydim::examples::Elliptic_PCC_3D::test::Test_Types TestType() const
    { return (Polydim::examples::Elliptic_PCC_3D::test::Test_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("TestType"); }
    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D MeshGenerator() const
    { return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"); }
    inline std::string MeshImportFilePath() const
    { return Gedim::Configurations::GetPropertyValue<string>("MeshImportFilePath"); }
    inline double MeshMaxVolume() const
    { return Gedim::Configurations::GetPropertyValue<double>("MeshMaxVolume"); }
    inline double GeometricTolerance1D() const
    { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D"); }
    inline double GeometricTolerance2D() const
    { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance2D"); }
    inline double GeometricTolerance3D() const
    { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance3D"); }

    inline Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Types VemType() const
    { return (Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("VemType"); }
    inline bool ComputeVEMPerformance() const
    { return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance"); }
    inline unsigned int VemOrder() const
    { return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder"); }

};
}
}
}

#endif
