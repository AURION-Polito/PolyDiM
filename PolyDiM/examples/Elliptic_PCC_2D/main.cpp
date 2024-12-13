#include "MeshUtilities.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"
#include "Eigen_CholeskySolver.hpp"
#include "assembler.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"

struct PatchTest final
{
    static unsigned int order;

    static Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain()
    {
      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

      domain.area = 1.0;

      domain.vertices = Eigen::MatrixXd::Zero(3, 4);
      domain.vertices.row(0)<< 0.0, 1.0, 1.0, 0.0;
      domain.vertices.row(1)<< 0.0, 0.0, 1.0, 1.0;

      return domain;
    }

    static std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info()
    {
      return {
        { 0, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
        { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 }  },
        { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } }
      };
    }

    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
    { return Eigen::VectorXd::Constant(points.cols(), 1.0); };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
      Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(),
                                                              2.0 * order * (order - 1));
      const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

      const int max_order = order - 2;
      for (int i = 0; i < max_order; ++i)
        source_term.array() *= polynomial;

      return - source_term;
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
      if (marker != 1)
        throw std::runtime_error("Unknown marker");

      const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

      Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
      for (int i = 0; i < order; ++i)
        result.array() *= polynomial;

      return result;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    { throw std::runtime_error("Not supported"); }

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {

      const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

      Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
      for (int i = 0; i < order; ++i)
        result.array() *= polynomial;

      return result;
    };

    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
    {
      Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
      const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

      const int max_order = order - 1;
      for (int i = 0; i < max_order; ++i)
        derivatives.array() *= polynomial;

      return
      {
        derivatives,
            derivatives,
            Eigen::VectorXd::Zero(points.cols())
      };
    }
};

struct Poisson_Polynomial_Problem final
{
    static Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain()
    {
      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

      domain.area = 1.0;

      domain.vertices = Eigen::MatrixXd::Zero(3, 4);
      domain.vertices.row(0)<< 0.0, 1.0, 1.0, 0.0;
      domain.vertices.row(1)<< 0.0, 0.0, 1.0, 1.0;

      return domain;
    }

    static std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info()
    {
      return {
        { 0, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
        { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 }  },
        { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } },
        { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
        { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4 } }
      };
    }

    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
    {
      const double k = 1.0;
      return Eigen::VectorXd::Constant(points.cols(), k);
    };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
      return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                     points.row(0).array() * (1.0 - points.row(0).array()));
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
      if (marker != 1)
        throw std::runtime_error("Unknown marker");

      return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                     points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    {
      switch(marker)
      {
        case 2: // co-normal derivatives on the right
          return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 4: // co-normal derivatives on the left
          return - 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
          throw std::runtime_error("Unknown marker");
      }
    }

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {
      return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                     points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    };

    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
    {
      return
      {
        16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
            16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
            Eigen::VectorXd::Zero(points.cols())
      };
    }
};

unsigned int PatchTest::order;

int main(int argc, char** argv)
{
  Elliptic_PCC_2D::Program_configuration config;

  if (!Gedim::Output::FileExists("./Parameters.ini"))
    Gedim::Configurations::ExportToIni("./Parameters.ini",
                                       false);
  else
    Gedim::Configurations::InitializeFromIni("./Parameters.ini");

  Gedim::Configurations::Initialize(argc, argv);

  Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
  geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance();
  Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

  Gedim::MeshUtilities meshUtilities;

  /// Create folders
  const string exportFolder = config.ExportFolder();
  Gedim::Output::CreateFolder(exportFolder);

  const string exportCsvFolder = exportFolder + "/Mesh";
  Gedim::Output::CreateFolder(exportCsvFolder);
  const string exportVtuFolder = exportFolder + "/Paraview";
  Gedim::Output::CreateFolder(exportVtuFolder);
  const string exportSolutionFolder = exportFolder + "/Solution";
  Gedim::Output::CreateFolder(exportSolutionFolder);

  const string logFolder = exportFolder + "/Log";

  /// Set Profiler
  Gedim::Profiler::ActivateProfiler = true;

  /// Set Log folder
  Gedim::Output::CreateFolder(logFolder);
  Gedim::LogFile::LogFolder = logFolder;

  /// Export Configuration of the following Run
  Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini",
                                     false);

  /// Create domain
  Gedim::Output::PrintGenericMessage("CreateDomain...", true);
  Gedim::Profiler::StartTime("CreateDomain");

  PatchTest::order = config.VemOrder();

  const auto domain = PatchTest::domain();
  const auto boundary_info = PatchTest::boundary_info();
  const auto diffusion_term = PatchTest::diffusion_term;
  const auto source_term = PatchTest::source_term;
  const auto weak_boundary_condition = PatchTest::weak_boundary_condition;
  const auto strong_boundary_condition = PatchTest::strong_boundary_condition;
  const auto exact_solution = PatchTest::exact_solution;
  const auto exact_derivative_solution = PatchTest::exact_derivative_solution;

  Gedim::Profiler::StopTime("CreateDomain");
  Gedim::Output::PrintStatusProgram("CreateDomain");

  // export domain
  {
    Gedim::VTKUtilities vtkUtilities;
    vtkUtilities.AddPolygon(domain.vertices);
    vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
  }

  /// Create domain mesh
  Gedim::Output::PrintGenericMessage("CreateMesh...", true);
  Gedim::Profiler::StartTime("CreateMesh");

  Gedim::MeshMatrices meshData;
  Gedim::MeshMatricesDAO mesh(meshData);

  switch (config.MeshGenerator())
  {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    {
      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                  meshUtilities,
                                                                  config.MeshGenerator(),
                                                                  domain,
                                                                  config.MeshMaxArea(),
                                                                  mesh);
    }
      break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter:
    {
      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(geometryUtilities,
                                                                  meshUtilities,
                                                                  config.MeshGenerator(),
                                                                  config.MeshImportFilePath(),
                                                                  mesh);
    }
      break;
    default:
      throw runtime_error("MeshGenerator " +
                          to_string((unsigned int)config.MeshGenerator()) +
                          " not supported");
  }

  Gedim::Profiler::StopTime("CreateMesh");
  Gedim::Output::PrintStatusProgram("CreateMesh");

  // Export the domain mesh
  {
    meshUtilities.ExportMeshToVTU(mesh,
                                  exportVtuFolder,
                                  "Domain_Mesh");
  }

  Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
  Gedim::Profiler::StartTime("ComputeGeometricProperties");

  const auto meshGeometricData = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities,
                                                                                                            meshUtilities,
                                                                                                            mesh);
  Gedim::Profiler::StopTime("ComputeGeometricProperties");
  Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

  /// Initialize Discrete Space
  Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
  Gedim::Profiler::StartTime("CreateVEMSpace");

  Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;
  const auto reference_element_data = vem_reference_element.Create(config.VemOrder());

  Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data =
  { mesh };

  Polydim::PDETools::DOFs::DOFsManager dofManager;
  const auto meshDOFsInfo = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                                   {
                                                                     {
                                                                       reference_element_data.NumDofs0D,
                                                                       reference_element_data.NumDofs1D,
                                                                       reference_element_data.NumDofs2D,
                                                                       0
                                                                     },
                                                                     boundary_info
                                                                   });

  const auto dofs_data = dofManager.CreateDOFs<2>(meshDOFsInfo,
                                                  mesh_connectivity_data);

  Gedim::Output::PrintGenericMessage("VEM Space with " +
                                     to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                     to_string(dofs_data.NumberStrongs) + " STRONGs", true);

  Gedim::Profiler::StopTime("CreateVEMSpace");
  Gedim::Output::PrintStatusProgram("CreateVEMSpace");

  Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
  Gedim::Profiler::StartTime("AssembleSystem");

  Elliptic_PCC_2D::Assembler assembler;

  auto assembler_data = assembler.Assemble(geometryUtilities,
                                           mesh,
                                           meshGeometricData,
                                           meshDOFsInfo,
                                           dofs_data,
                                           reference_element_data,
                                           diffusion_term,
                                           source_term,
                                           strong_boundary_condition,
                                           weak_boundary_condition);

  Gedim::Profiler::StopTime("AssembleSystem");
  Gedim::Output::PrintStatusProgram("AssembleSystem");

  if (dofs_data.NumberDOFs > 0)
  {
    Gedim::Output::PrintGenericMessage("Factorize...", true);
    Gedim::Profiler::StartTime("Factorize");


    Gedim::Eigen_CholeskySolver solver;
    solver.Initialize(assembler_data.globalMatrixA);

    Gedim::Profiler::StopTime("Factorize");
    Gedim::Output::PrintStatusProgram("Factorize");

    Gedim::Output::PrintGenericMessage("Solve...", true);
    Gedim::Profiler::StartTime("Solve");

    solver.Solve(assembler_data.rightHandSide,
                 assembler_data.solution);

    Gedim::Profiler::StopTime("Solve");
    Gedim::Output::PrintStatusProgram("Solve");
  }

  Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
  Gedim::Profiler::StartTime("ComputeErrors");

  auto post_process_data = assembler.PostProcessSolution(geometryUtilities,
                                                         mesh,
                                                         meshGeometricData,
                                                         dofs_data,
                                                         reference_element_data,
                                                         assembler_data,
                                                         exact_solution,
                                                         exact_derivative_solution);

  Gedim::Profiler::StopTime("ComputeErrors");
  Gedim::Output::PrintStatusProgram("ComputeErrors");

  Gedim::Output::PrintGenericMessage("ExportSolution...", true);
  Gedim::Profiler::StartTime("ExportSolution");

  {
    const char separator = ';';

    std::cout<< "VemType" << separator;
    std::cout<< "VemOrder" << separator;
    std::cout<< "Cell2Ds" <<  separator;
    std::cout<< "Dofs" <<  separator;
    std::cout<< "Strongs" <<  separator;
    std::cout<< "h" <<  separator;
    std::cout<< "errorL2" <<  separator;
    std::cout<< "errorH1" << separator;
    std::cout<< "normL2" <<  separator;
    std::cout<< "normH1" << separator;
    std::cout<< "nnzA" << separator;
    std::cout<< "residual" << std::endl;

    std::cout.precision(2);
    std::cout<< scientific<< static_cast<unsigned int>(config.VemType())<< separator;
    std::cout<< scientific<< config.VemOrder()<< separator;
    std::cout<< scientific<< mesh.Cell2DTotalNumber()<< separator;
    std::cout<< scientific<< dofs_data.NumberDOFs<< separator;
    std::cout<< scientific<< dofs_data.NumberStrongs<< separator;
    std::cout<< scientific<< post_process_data.mesh_size << separator;
    std::cout<< scientific<< post_process_data.error_L2<< separator;
    std::cout<< scientific<< post_process_data.error_H1<< separator;
    std::cout<< scientific<< post_process_data.norm_L2<< separator;
    std::cout<< scientific<< post_process_data.norm_H1<< separator;
    std::cout<< scientific<< assembler_data.globalMatrixA.NonZeros()<< separator;
    std::cout<< scientific<< post_process_data.residual_norm<< std::endl;
  }

  {
    const char separator = ';';
    const string errorFileName = exportSolutionFolder +
                                 "/Errors.csv";
    const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

    std::ofstream errorFile(errorFileName,
                            std::ios_base::app | std::ios_base::out);
    if (!errorFileExists)
    {
      errorFile<< "VemType" << separator;
      errorFile<< "VemOrder" << separator;
      errorFile<< "Cell2Ds" <<  separator;
      errorFile<< "Dofs" <<  separator;
      errorFile<< "Strongs" <<  separator;
      errorFile<< "h" <<  separator;
      errorFile<< "errorL2" <<  separator;
      errorFile<< "errorH1" << separator;
      errorFile<< "normL2" <<  separator;
      errorFile<< "normH1" << separator;
      errorFile<< "nnzA" << separator;
      errorFile<< "residual" << std::endl;
    }

    errorFile.precision(16);
    errorFile<< scientific<< static_cast<unsigned int>(config.VemType())<< separator;
    errorFile<< scientific<< config.VemOrder()<< separator;
    errorFile<< scientific<< mesh.Cell2DTotalNumber()<< separator;
    errorFile<< scientific<< dofs_data.NumberDOFs<< separator;
    errorFile<< scientific<< dofs_data.NumberStrongs<< separator;
    errorFile<< scientific<< post_process_data.mesh_size << separator;
    errorFile<< scientific<< post_process_data.error_L2<< separator;
    errorFile<< scientific<< post_process_data.error_H1<< separator;
    errorFile<< scientific<< post_process_data.norm_L2<< separator;
    errorFile<< scientific<< post_process_data.norm_H1<< separator;
    errorFile<< scientific<< assembler_data.globalMatrixA.NonZeros()<< separator;
    errorFile<< scientific<< post_process_data.residual_norm<< std::endl;

    errorFile.close();
  }

  {
    {
      Gedim::VTKUtilities exporter;
      exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                           mesh.Cell2DsVertices(),
                           {
                             {
                               "Numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
                               post_process_data.cell0Ds_numeric.data()
                             },
                             {
                               "Exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
                               post_process_data.cell0Ds_exact.data()
                             },
                             {
                               "ErrorL2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
                               post_process_data.cell2Ds_error_L2.data()
                             },
                             {
                               "ErrorH1",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
                               post_process_data.cell2Ds_error_H1.data()
                             }
                           });

      exporter.Export(exportVtuFolder + "/Solution_Cell2Ds.vtu");
    }
  }

  Gedim::Profiler::StopTime("ExportSolution");
  Gedim::Output::PrintStatusProgram("ExportSolution");

  Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
  Gedim::Profiler::StartTime("ComputeVEMPerformance");

  if (config.ComputeVEMPerformance())
  {
    const auto vemPerformance = assembler.ComputeVemPerformance(geometryUtilities,
                                                                mesh,
                                                                meshGeometricData,
                                                                reference_element_data);
    {
      const char separator = ',';
      /// Export Cell2Ds VEM performance
      ofstream exporter;

      exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance.csv");
      exporter.precision(16);

      if (exporter.fail())
        throw runtime_error("Error on mesh cell2Ds file");

      exporter<< "Cell2D_Index"<< separator;
      exporter<< "NumQuadPoints_Boundary"<< separator;
      exporter<< "NumQuadPoints_Internal"<< separator;
      exporter<< "PiNabla_Cond"<< separator;
      exporter<< "Pi0k_Cond"<< separator;
      exporter<< "Pi0km1_Cond"<< separator;
      exporter<< "PiNabla_Error"<< separator;
      exporter<< "Pi0k_Error"<< separator;
      exporter<< "Pi0km1_Error"<< separator;
      exporter<< "HCD_Error"<< separator;
      exporter<< "GBD_Error"<< separator;
      exporter<< "Stab_Error"<< endl;

      for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
      {
        const auto& cell2DPerformance = vemPerformance.Cell2DsPerformance[v].Analysis;

        exporter<< scientific<< v<< separator;
        exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumBoundaryQuadraturePoints<< separator;
        exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumInternalQuadraturePoints<< separator;
        exporter<< scientific<< cell2DPerformance.PiNablaConditioning<< separator;
        exporter<< scientific<< cell2DPerformance.Pi0kConditioning<< separator;
        exporter<< scientific<< cell2DPerformance.Pi0km1Conditioning<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorPiNabla<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorPi0k<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorPi0km1<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorHCD<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorGBD<< separator;
        exporter<< scientific<< cell2DPerformance.ErrorStabilization<< endl;
      }

      exporter.close();
    }
  }

  Gedim::Profiler::StopTime("ComputeVEMPerformance");
  Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

  return 0;
}
