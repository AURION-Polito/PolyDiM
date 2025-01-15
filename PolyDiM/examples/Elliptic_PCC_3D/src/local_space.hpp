#ifndef __local_space_H
#define __local_space_H

#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"
#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_3D_Creator.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
namespace local_space
{
struct ReferenceElement_Data final
{
    Program_configuration::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement> VEM_ReferenceElement_2D;
    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement_Data_2D;
    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_3D_ReferenceElement> VEM_ReferenceElement_3D;
    Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data VEM_ReferenceElement_Data_3D;
    VEM::PCC::VEM_PCC_3D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::PCC::I_VEM_PCC_3D_LocalSpace> VEM_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement> FEM_ReferenceElement_3D;
    Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement_Data FEM_ReferenceElement_Data_3D;
    std::unique_ptr<FEM::PCC::FEM_Tetrahedron_PCC_3D_LocalSpace> FEM_LocalSpace;
};

struct LocalSpace_Data final
{
    struct VEM_Geometry final
    {
        std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> PolygonalFaces;
        Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry Polyhedron;
    };

    VEM_Geometry VEM_Geometry;
    Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Data VEM_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_LocalSpace_Data FEM_LocalSpace_Data;
};

struct Performance_Data final
{
    struct Cell3D_Performance final
    {
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
    };

    Cell3D_Performance VEM_Performance_Data;
};

ReferenceElement_Data CreateReferenceElement(const Program_configuration::MethodTypes &method_type, const unsigned int method_order);

std::array<unsigned int, 4> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data);

LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                 const unsigned int cell3D_index,
                                 const ReferenceElement_Data &reference_element_data);

Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                     const LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

Eigen::MatrixXd BasisFunctionsValuesOnFace(const unsigned int &face_local_index,
                                           const ReferenceElement_Data &reference_element_data,
                                           const LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &quadrature_points);

std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

Eigen::MatrixXd StabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::PiNabla);

Eigen::MatrixXd EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index);

Gedim::Quadrature::QuadratureData FaceDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                                      const LocalSpace_Data &local_space_data,
                                                      const unsigned int face_local_index,
                                                      unsigned int &quadrature_offset);

Eigen::VectorXd FaceDofs(const ReferenceElement_Data &reference_element_data,
                         const LocalSpace_Data &local_space_data,
                         const unsigned int face_local_index,
                         const Eigen::VectorXd &strong_values,
                         const Gedim::Quadrature::QuadratureData &quadrature_offset);

Gedim::Quadrature::QuadratureData FaceQuadrature(const ReferenceElement_Data &reference_element_data,
                                                 const LocalSpace_Data &local_space_data,
                                                 const unsigned int face_local_index,
                                                 unsigned int &quadrature_offset);

Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data);

unsigned int Size(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

}; // namespace local_space
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
