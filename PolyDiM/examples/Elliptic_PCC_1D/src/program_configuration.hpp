#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "FEM_PCC_1D_Creator.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_1D
{
struct Program_configuration final
{
  public:
    Program_configuration()
    {
        Gedim::Configurations::AddProperty("TestType",
                                           static_cast<unsigned int>(Polydim::examples::Elliptic_PCC_1D::test::Test_Types::Patch_Test),
                                           "Test Type 1 - Patch_Test; 2 - Poisson_Polynomial_Problem "
                                           "(Default: 1)");
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Run", "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty(
            "MeshGenerator",
            static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Equispaced),
            "Mesh 1D gereator type, 0 - Equispaced; 1 - Imported; 2 - Minimal (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxLength", 0.1, "Mesh 1D maximum relative cell length (Default: 0.1)");

        Gedim::Configurations::AddProperty("GeometricTolerance1D", 1.0e-12, "Geometric Tolerance 1D (Default: 1.0e-12)");

        // Method parameters
        Gedim::Configurations::AddProperty("MethodType",
                                           static_cast<unsigned int>(Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace),
                                           "Method Type, 1 - FEM (Default: 1)");
        Gedim::Configurations::AddProperty("MethodOrder", static_cast<unsigned int>(1), "Method order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeMethodPerformance", true, "Compute Method Performance (Default: true)");
    }

    inline string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<string>("ExportFolder");
    }

    inline Polydim::examples::Elliptic_PCC_1D::test::Test_Types TestType() const
    {
        return (Polydim::examples::Elliptic_PCC_1D::test::Test_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("TestType");
    }
    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D MeshGenerator() const
    {
        return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D)
            Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator");
    }
    inline std::string MeshImportFilePath() const
    {
        return Gedim::Configurations::GetPropertyValue<string>("MeshImportFilePath");
    }
    inline double MeshMaxLength() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("MeshMaxLength");
    }
    inline double GeometricTolerance1D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D");
    }

    inline bool ComputeMethodPerformance() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeMethodPerformance");
    }
    inline unsigned int MethodOrder() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("MethodOrder");
    }
    inline Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Types MethodType() const
    {
        return (Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("M"
                                                                                                                     "e"
                                                                                                                     "t"
                                                                                                                     "h"
                                                                                                                     "o"
                                                                                                                     "d"
                                                                                                                     "T"
                                                                                                                     "y"
                                                                                                                     "p"
                                                                                                                     "e");
    }
};
} // namespace Elliptic_PCC_1D
} // namespace examples
} // namespace Polydim

#endif
