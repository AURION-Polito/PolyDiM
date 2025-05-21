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

#ifndef __local_space_H
#define __local_space_H

#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_2D_Creator.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
namespace examples
{
namespace Elastic_PCC_2D
{
namespace local_space
{
struct ReferenceElement_Data final
{
    Program_configuration::MethodTypes Method_Type;
    unsigned int Order;
    unsigned int Dimension;

    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement> VEM_ReferenceElement;
    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement_Data;
    VEM::PCC::VEM_PCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::PCC::I_VEM_PCC_2D_LocalSpace> VEM_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement> FEM_ReferenceElement;
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
    std::unique_ptr<FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace> FEM_LocalSpace;
};

struct LocalSpace_Data final
{
    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data VEM_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_Polygon_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data FEM_LocalSpace_Data;
};

struct Performance_Data final
{
    struct Cell2D_Performance final
    {
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
    };

    Cell2D_Performance VEM_Performance_Data;
};

ReferenceElement_Data CreateReferenceElement(const Program_configuration::MethodTypes &method_type, const unsigned int method_order);

std::array<unsigned int, 4> ReferenceElementNumDOFs(const ReferenceElement_Data &reference_element_data);

LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data);

std::vector<Eigen::MatrixXd> BasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0k);

Eigen::MatrixXd BasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                            const ReferenceElement_Data &reference_element_data,
                                            const LocalSpace_Data &local_space_data,
                                            const Eigen::MatrixXd &pointsCurvilinearCoordinates);

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

Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data);

unsigned int Size(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data);

}; // namespace local_space
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
