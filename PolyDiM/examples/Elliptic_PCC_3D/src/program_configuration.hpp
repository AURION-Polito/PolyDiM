// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
struct Program_configuration final
{
    enum struct MethodTypes
    {
        FEM_PCC = 0,
        VEM_PCC = 1,
        VEM_PCC_Inertia = 2,
        VEM_PCC_Ortho = 3
    };

  public:
    Program_configuration()
    {
        Gedim::Configurations::AddProperty("TestType",
                                           static_cast<unsigned int>(Polydim::examples::Elliptic_PCC_3D::test::Test_Types::Patch_Test),
                                           "Test Type 1 - Patch_Test; 2 - Poisson_Polynomial_Problem "
                                           "(Default: 1)");
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Run", "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty(
            "MeshGenerator",
            static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Tetrahedral),
            "Mesh 3D gereator type, 0 - Tetrahedral; 1 - Minimal; 2 - "
            "Polyhedral; 3 - OVMImporter; 4 - VtkImporter; 5 - CsvImporter; 6 - Cubic "
            "(Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxVolume", 0.1, "Mesh 3D maximum relative cell volume (Default: 0.1)");

        Gedim::Configurations::AddProperty("GeometricTolerance1D", 1.0e-12, "Geometric Tolerance 1D (Default: 1.0e-12)");

        Gedim::Configurations::AddProperty("GeometricTolerance2D", 1.0e-14, "Geometric Tolerance 2D (Default: 1.0e-14)");

        Gedim::Configurations::AddProperty("GeometricTolerance3D", 1.0e-15, "Geometric Tolerance 3D (Default: 1.0e-15)");

        // Method parameters
        Gedim::Configurations::AddProperty("MethodType",
                                           static_cast<unsigned int>(MethodTypes::FEM_PCC),
                                           "Method Type, 0 - FEM_PCC; 1 - EVem; 2 - EVem_Inertia; 3 - "
                                           "EVem_Ortho (Default: "
                                           "0)");
        Gedim::Configurations::AddProperty("MethodOrder", static_cast<unsigned int>(1), "Method order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeMethodPerformance", true, "Compute Method Performance (Default: true)");
    }

    inline std::string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("ExportFolder");
    }

    inline Polydim::examples::Elliptic_PCC_3D::test::Test_Types TestType() const
    {
        return (Polydim::examples::Elliptic_PCC_3D::test::Test_Types)Gedim::Configurations::GetPropertyValue<unsigned int>("TestType");
    }
    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D MeshGenerator() const
    {
        return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D)
            Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator");
    }
    inline std::string MeshImportFilePath() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("MeshImportFilePath");
    }
    inline double MeshMaxVolume() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("MeshMaxVolume");
    }
    inline double GeometricTolerance1D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D");
    }
    inline double GeometricTolerance2D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance2D");
    }
    inline double GeometricTolerance3D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance3D");
    }

    inline MethodTypes MethodType() const
    {
        return (MethodTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("MethodType");
    }
    inline bool ComputeMethodPerformance() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeMethodPerformance");
    }
    inline unsigned int MethodOrder() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("MethodOrder");
    }
};
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
