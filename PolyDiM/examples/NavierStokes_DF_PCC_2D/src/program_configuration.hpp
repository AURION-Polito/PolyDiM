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
namespace NavierStokes_DF_PCC_2D
{
struct Program_configuration final
{
    enum ConvectiveFormType
    {
        Conv = 0,
        Skew = 1,
        None = 2
    };

    Program_configuration()
    {
        Gedim::Configurations::AddProperty(
            "TestType",
            static_cast<unsigned int>(Polydim::examples::NavierStokes_DF_PCC_2D::test::Test_Types::Patch_Test),
            "Test Type 1 - Patch_Test; 2 - StokesSinSin; 3 - NavierStokes; 4 - NavierStokes_VanishingExternalLoad "
            "(Default: 1)");

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
            "Vem Type, 1 - Vem; (Default: 1)");
        Gedim::Configurations::AddProperty("VemOrder", static_cast<unsigned int>(2), "VEM order (Default: 2)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance", true, "Compute VEM Performance (Default: true)");

        /// Solver
        Gedim::Configurations::AddProperty("ConvectiveForm",
                                           static_cast<unsigned int>(ConvectiveFormType::Conv),
                                           "Convective Form, 0 - Conv; 1 - Skew; 2 - None (Default: 0)");

        Gedim::Configurations::AddProperty("NLMaxNumberIterations",
                                           static_cast<unsigned int>(1),
                                           "Maximum number of non-linear iterations (Default: 1)");

        Gedim::Configurations::AddProperty("NLAbsChangeInSolutionTolerance",
                                           1.0e-06,
                                           "Absolute tolerance for the non-linear method - change in solution "
                                           "(Default: 1.0e-06)");

        Gedim::Configurations::AddProperty("NLAbsResidualTolerance",
                                           1.0e-06,
                                           "Absolute tolerance for the non-linear method - residual (Default: "
                                           "1.0e-06)");

        Gedim::Configurations::AddProperty("NLRelChangeInSolutionTolerance",
                                           1.0e-06,
                                           "Relative tolerance for the non-linear method - change in solution "
                                           "(Default: 1.0e-06)");

        Gedim::Configurations::AddProperty("NLRelResidualTolerance",
                                           1.0e-06,
                                           "Relative tolerance for the non-linear method - residual (Default: "
                                           "1.0e-06)");
    }

    inline Polydim::examples::NavierStokes_DF_PCC_2D::test::Test_Types TestType() const
    {
        return (Polydim::examples::NavierStokes_DF_PCC_2D::test::Test_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("TestType");
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
        if (Gedim::Configurations::GetPropertyValue<unsigned int>("VemType") == 2)
            throw runtime_error("not valid vem type");

        return (Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("VemType");
    }
    inline bool ComputeVEMPerformance() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance");
    }
    inline unsigned int VemOrder() const
    {
        if (Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder") < 2)
            throw runtime_error("not valid order");
        return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder");
    }

    inline ConvectiveFormType ConvectiveForm() const
    {
        return (ConvectiveFormType)Gedim::Configurations::GetPropertyValue<unsigned int>("ConvectiveForm");
    }
    inline unsigned int NLMaxNumberIterations() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("NLMaxNumberIterations");
    }
    inline double NLAbsChangeInSolutionTolerance() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("NLAbsChangeInSolutionTolerance");
    }
    inline double NLAbsResidualTolerance() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("NLAbsResidualTolerance");
    }
    inline double NLRelChangeInSolutionTolerance() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("NLRelChangeInSolutionTolerance");
    }
    inline double NLRelResidualTolerance() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("NLRelResidualTolerance");
    }
};
} // namespace NavierStokes_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
