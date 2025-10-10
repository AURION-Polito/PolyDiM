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

#include "program_utilities.hpp"

#include "VTKUtilities.hpp"
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Parabolic_PCC_2D::test::I_Test> create_test(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Patch_Test>();
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Elliptic_Polynomial_Problem:
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Elliptic_Polynomial_Problem>();
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Elliptic_Problem:
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Elliptic_Problem>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::RandomDistorted: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(meshUtilities, config.MeshGenerator(), config.MeshImportFilePath(), mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator " + std::to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities, meshUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Static_Problem_Data &assembler_static_data,
                     const Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';

        std::cout << "ProgramType" << separator;
        std::cout << "MethodType" << separator;
        std::cout << "MethodOrder" << separator;
        std::cout << "Cell2Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorL2" << separator;
        std::cout << "errorH1" << separator;
        std::cout << "normL2" << separator;
        std::cout << "normH1" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << Method_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << dofs_data.NumberDOFs << separator;
        std::cout << std::scientific << dofs_data.NumberStrongs << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_L2 << separator;
        std::cout << std::scientific << post_process_data.error_H1 << separator;
        std::cout << std::scientific << post_process_data.norm_L2 << separator;
        std::cout << std::scientific << post_process_data.norm_H1 << separator;
        std::cout << std::scientific << assembler_static_data.globalMatrixA.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(Method_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "MethodOrder" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorL2" << separator;
            errorFile << "errorH1" << separator;
            errorFile << "normL2" << separator;
            errorFile << "normH1" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << Method_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << dofs_data.NumberDOFs << separator;
        errorFile << std::scientific << dofs_data.NumberStrongs << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_L2 << separator;
        errorFile << std::scientific << post_process_data.error_H1 << separator;
        errorFile << std::scientific << post_process_data.norm_L2 << separator;
        errorFile << std::scientific << post_process_data.norm_H1 << separator;
        errorFile << std::scientific << assembler_static_data.globalMatrixA.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                                 mesh.Cell2DsVertices(),
                                 {{"Numeric",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
                                   post_process_data.cell0Ds_numeric.data()},
                                  {"Exact",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
                                   post_process_data.cell0Ds_exact.data()},
                                  {"ErrorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
                                   post_process_data.cell2Ds_error_L2.data()},
                                  {"ErrorH1",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
                                   post_process_data.cell2Ds_error_H1.data()}});

            exporter.Export(exportVtuFolder + "/Solution_" + std::to_string(TEST_ID) + "_" + std::to_string(Method_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_dofs(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Problem_Data &assembler_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                 const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> solution_values;
    std::list<double> rhs_values;
    std::list<double> dof_global_index_values;
    std::list<double> dof_type_values;
    std::list<double> dof_cell_index_values;
    std::list<double> dof_dimension_values;
    std::list<double> dof_boundary_type_values;
    std::list<double> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[0].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(0);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[1].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const std::vector<double> local_edge_coordinates = geometryUtilities.EquispaceCoordinates(num_loc_dofs, 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[2].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const auto local_polygon_coordinates = geometryUtilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell2DsCentroids.at(c);
        const auto polygonCentroidEdgesDistance =
            geometryUtilities.PolygonCentroidEdgesDistance(mesh_geometric_data.Cell2DsVertices.at(c),
                                                           mesh_geometric_data.Cell2DsCentroids.at(c),
                                                           mesh_geometric_data.Cell2DsEdgeNormals.at(c));
        const double circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);

            dofs_coordinate.push_back(
                polygon_centroid +
                circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  0.0));

            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    {
        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
        unsigned int c = 0;
        for (const auto &dof_coordinate : dofs_coordinate)
            coordinates.col(c++) << dof_coordinate;
        const auto rhs_values_data = std::vector<double>(rhs_values.begin(), rhs_values.end());
        const auto solution_values_data = std::vector<double>(solution_values.begin(), solution_values.end());
        const auto dof_global_index_values_data =
            std::vector<double>(dof_global_index_values.begin(), dof_global_index_values.end());
        const auto dof_type_values_data = std::vector<double>(dof_type_values.begin(), dof_type_values.end());
        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
        const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
        const auto dof_boundary_type_values_data =
            std::vector<double>(dof_boundary_type_values.begin(), dof_boundary_type_values.end());
        const auto dof_boundary_marker_values_data =
            std::vector<double>(dof_boundary_marker_values.begin(), dof_boundary_marker_values.end());

        Gedim::VTKUtilities exporter;
        exporter.AddPoints(coordinates,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values_data.size()),
                             dof_dimension_values_data.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values_data.size()),
                             dof_cell_index_values_data.data()},
                            {"boundary_type",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
                             dof_boundary_type_values_data.data()},
                            {"boundary_marker",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
                             dof_boundary_marker_values_data.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()},
                            {"dof_type",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_type_values_data.size()),
                             dof_type_values_data.data()},
                            {"rhs",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(rhs_values_data.size()),
                             rhs_values_data.data()},
                            {"solution",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(solution_values_data.size()),
                             solution_values_data.data()}});

        const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.Export(exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" + std::to_string(Method_ID) + +"_" +
                        std::to_string(config.MethodOrder()) + ".vtu");
    }
}
// ***************************************************************************
void export_performance(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                        const Assembler::Performance_Data &performance_data,
                        const std::string &exportFolder)
{
    {
        const char separator = ',';
        std::ofstream exporter;
        const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.open(exportFolder + "/Cell2Ds_MethodPerformance_" + std::to_string(TEST_ID) + "_" +
                      std::to_string(Method_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv");
        exporter.precision(16);

        if (exporter.fail())
            throw std::runtime_error("Error on mesh cell2Ds file");

        exporter << "Cell2D_Index" << separator;
        exporter << "NumQuadPoints_Boundary" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "PiNabla_Cond" << separator;
        exporter << "Pi0k_Cond" << separator;
        exporter << "Pi0km1_Cond" << separator;
        exporter << "PiNabla_Error" << separator;
        exporter << "Pi0k_Error" << separator;
        exporter << "Pi0km1_Error" << separator;
        exporter << "HCD_Error" << separator;
        exporter << "GBD_Error" << separator;
        exporter << "Stab_Error" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell2DsPerformance.size(); v++)
        {
            const auto &cell2D_performance = performance_data.Cell2DsPerformance[v].VEM_Performance_Data;

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell2D_performance.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.Analysis.PiNablaConditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.Pi0kConditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.Pi0km1Conditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPiNabla << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPi0k << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPi0km1 << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorHCD << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorGBD << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorStabilization << std::endl;
        }

        exporter.close();
    }
}
// ***************************************************************************
std::vector<double> create_time_steps(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                      const std::array<double, 2>& time_domain)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::Output::Assert(geometryUtilities.IsValuePositive(config.TimeStep(), geometryUtilities.Tolerance1D()));
    std::vector<double> times = geometryUtilities.EquispaceCoordinates(config.TimeStep(), true);
    for (unsigned int t = 0; t < times.size(); t++)
        times[t] = (time_domain.at(1) - time_domain.at(0)) * times[t] + time_domain.at(0);

    return times;
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim
