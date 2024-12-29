#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "VEM_DF_PCC_2D_Creator.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
struct Program_configuration final
{
    Program_configuration()
    {
        Gedim::Configurations::AddProperty("TestType",
                                           static_cast<unsigned int>(Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Patch_Test),
                                           "Test Type: 1 - Patch_Test; 2 - StokesSinSin; 3 - Stokes_ZeroVelocity_1; 4 "
                                           "- Stokes_ZeroVelocity_2; 5 - Darcy; 6 - Brinkman (Default: 1)");

        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Run", "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty(
            "MeshGenerator",
            static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular),
            "Mesh 2D gereator type, 0 - Triangular; 1 - Minimal; 2 - "
            "Polygonal; 3 - OFF Importer; 4 - Csv Importer (semicolon); 5 - Squared (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxArea", 0.1, "Mesh 2D maximum relative cell area (Default: 0.1)");

        Gedim::Configurations::AddProperty("GeometricTolerance1D", 1.0e-12, "Geometric Tolerance 1D (Default: 1.0e-12)");

        Gedim::Configurations::AddProperty("GeometricTolerance2D", 1.0e-14, "Geometric Tolerance 2D (Default: 1.0e-14)");

        /// Method parameters
        Gedim::Configurations::AddProperty(
            "VemType",
            static_cast<unsigned int>(Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types::VEM_DF_PCC_2D_LocalSpace),
            "Vem Type, 1 - Vem; 2 - RVem (Default: 1)");
        Gedim::Configurations::AddProperty("VemOrder", static_cast<unsigned int>(2), "VEM order (Default: 2)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance", true, "Compute VEM Performance (Default: true)");
        Gedim::Configurations::AddProperty("ComputeDiscrepancyError", false, "Compute Discrepancy error (Default: false)");
    }

    inline Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types TestType() const
    {
        return (Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("TestType");
    }

    inline string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<string>("ExportFolder");
    }

    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D MeshGenerator() const
    {
        return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D)
            Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator");
    }
    inline std::string MeshImportFilePath() const
    {
        return Gedim::Configurations::GetPropertyValue<string>("MeshImportFilePath");
    }
    inline double MeshMaxArea() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("MeshMaxArea");
    }
    inline double GeometricTolerance1D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D");
    }
    inline double GeometricTolerance2D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance2D");
    }

    inline Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types VemType() const
    {
        return (Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("VemType");
    }
    inline bool ComputeVEMPerformance() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance");
    }
    inline bool ComputeDiscrepancyError() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeDiscrepancyError");
    }
    inline unsigned int VemOrder() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder");
    }
};
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
