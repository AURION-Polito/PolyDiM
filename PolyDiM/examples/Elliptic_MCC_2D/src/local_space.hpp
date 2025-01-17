#ifndef __local_space_H
#define __local_space_H

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "VEM_MCC_2D_Creator.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{
namespace local_space
{
struct ReferenceElement_Data final
{
    Program_configuration::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Velocity_ReferenceElement> VEM_ReferenceElement_Velocity;
    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Pressure_ReferenceElement> VEM_ReferenceElement_Pressure;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data VEM_ReferenceElement_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data VEM_ReferenceElement_Data_Pressure;
    VEM::MCC::VEM_MCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Velocity_LocalSpace> VEM_LocalSpace_Velocity;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Pressure_LocalSpace> VEM_LocalSpace_Pressure;
};

struct LocalSpace_Data final
{
    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data VEM_LocalSpace_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data VEM_LocalSpace_Data_Pressure;
};

struct Performance_Data final
{
    struct Cell2D_Performance final
    {
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data Analysis;
    };

    Cell2D_Performance VEM_Performance_Data;
};

ReferenceElement_Data CreateReferenceElement(const Program_configuration::MethodTypes &method_type, const unsigned int method_order);

std::array<std::array<unsigned int, 4>, 2> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data);

LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const ReferenceElement_Data &reference_element_data,
                                                       const LocalSpace_Data &local_space_data);

Eigen::MatrixXd StabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

Gedim::Quadrature::QuadratureData EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                                      const LocalSpace_Data &local_space_data,
                                                      const unsigned int edge_local_index);

Eigen::VectorXd EdgeDofs(const ReferenceElement_Data &reference_element_data,
                         const LocalSpace_Data &local_space_data,
                         const unsigned int edge_local_index,
                         const Gedim::Quadrature::QuadratureData &edge_dofs_coordinates,
                         const Eigen::VectorXd &strong_values);

Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                                    const ReferenceElement_Data &reference_element_data,
                                                    const LocalSpace_Data &local_space_data,
                                                    const Eigen::MatrixXd &edge_quadrature_points);

Gedim::Quadrature::QuadratureData EdgeQuadrature(const ReferenceElement_Data &reference_element_data,
                                                 const LocalSpace_Data &local_space_data,
                                                 const unsigned int edge_local_index);

Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data);

unsigned int VelocitySize(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

}; // namespace local_space
} // namespace Elliptic_MCC_2D
} // namespace examples
} // namespace Polydim

#endif
