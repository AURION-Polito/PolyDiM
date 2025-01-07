#include "local_space.hpp"
#include <memory>


namespace Polydim
{
  namespace examples
  {
    namespace Elliptic_PCC_2D
    {
      namespace local_space
      {
        //***************************************************************************
        LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                         const unsigned int cell2D_index,
                                         const ReferenceElement_Data& reference_element_data)
        {
          LocalSpace_Data local_space_data;
          local_space_data.Method_Type = config.MethodType();

          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              local_space_data.FEM_Geometry = {
                config.GeometricTolerance1D(),
                config.GeometricTolerance2D(),
                mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index)
              };

              local_space_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace>();

              local_space_data.FEM_LocalSpace_Data = local_space_data.FEM_LocalSpace->CreateLocalSpace(reference_element_data.FEM_ReferenceElement,
                                                                                                       local_space_data.FEM_Geometry);
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              local_space_data.VEM_Geometry = {config.GeometricTolerance1D(),
                                               config.GeometricTolerance2D(),
                                               mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsCentroids.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsAreas.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsDiameters.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsTriangulations.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                               mesh_geometric_data.Cell2DsEdgeNormals.at(cell2D_index)};
              switch (local_space_data.Method_Type)
              {
                case Program_configuration::MethodTypes::VEM_PCC:
                  local_space_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace;
                  break;
                case Program_configuration::MethodTypes::VEM_PCC_Inertia:
                  local_space_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace;
                  break;
                case Program_configuration::MethodTypes::VEM_PCC_Ortho:
                  local_space_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Ortho_LocalSpace;
                  break;
                default:
                  throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
              }

              local_space_data.VEM_LocalSpace = Polydim::VEM::PCC::create_VEM_PCC_2D_local_space(local_space_data.VEM_Type);

              local_space_data.VEM_LocalSpace_Data = local_space_data.VEM_LocalSpace->CreateLocalSpace(reference_element_data.VEM_ReferenceElement,
                                                                                                       local_space_data.VEM_Geometry);
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }

          return local_space_data;
        }
        //***************************************************************************
        Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data& reference_element_data,
                                             const LocalSpace_Data& local_space_data)
        {
          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              return local_space_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement,
                                                                                  local_space_data.FEM_LocalSpace_Data);
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              return local_space_data.VEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data,
                                                                                  Polydim::VEM::PCC::ProjectionTypes::Pi0km1);
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }
        }
        //***************************************************************************
        std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(const ReferenceElement_Data& reference_element_data,
                                                                    const LocalSpace_Data& local_space_data)
        {
          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              return local_space_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement,
                                                                                            local_space_data.FEM_LocalSpace_Data);
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              return local_space_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_LocalSpace_Data,
                                                                                            Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }
        }
        //***************************************************************************
        Gedim::Quadrature::QuadratureData InternalQuadrature(const LocalSpace_Data& local_space_data)
        {
          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              return local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              return local_space_data.VEM_LocalSpace_Data.InternalQuadrature;
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }
        }
        //***************************************************************************
        unsigned int Size(const LocalSpace_Data& local_space_data)
        {
          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              return local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions;
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              return local_space_data.VEM_LocalSpace_Data.NumBasisFunctions;
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }
        }
        //***************************************************************************
        Eigen::MatrixXd StabilizationMatrix(const LocalSpace_Data& local_space_data)
        {
          switch (local_space_data.Method_Type)
          {
            case Program_configuration::MethodTypes::FEM_Triangle_PCC:
            {
              return Eigen::MatrixXd::Zero(local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                           local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
            }
            case Program_configuration::MethodTypes::VEM_PCC:
            case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            case Program_configuration::MethodTypes::VEM_PCC_Ortho:
            {
              return local_space_data.VEM_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data,
                                                                                         VEM::PCC::ProjectionTypes::PiNabla);
            }
            default:
              throw std::runtime_error("method type " + std::to_string((unsigned int)local_space_data.Method_Type) + " not supported");
          }
        }
        //***************************************************************************
     }
  } // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim
