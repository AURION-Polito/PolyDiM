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

#ifndef __FEM_PCC_2D_ReferenceElement_HPP
#define __FEM_PCC_2D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "FEM_Quadrilateral_PCC_2D_ReferenceElement.hpp"
#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "MeshUtilities.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data triangle_reference_element_data;
    Polydim::FEM::PCC::FEM_Quadrilateral_PCC_2D_ReferenceElement_Data quadrilateral_reference_element_data;

    Gedim::MeshUtilities::MeshGeometricData2DConfig mesh_geometric_data_config;

    FEM_PCC_2D_ReferenceElement_Data()
    {
        mesh_geometric_data_config.Cell2DsBoundingBox = false;
        mesh_geometric_data_config.Cell2DsTriangulations = false; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsAreas = false;          ///< cell2D areas
        mesh_geometric_data_config.Cell2DsCentroids = false;      ///< cell2D centroids
        mesh_geometric_data_config.Cell2DsDiameters = true;      ///< cell2D diameters
        mesh_geometric_data_config.Cell2DsEdgeDirections = false; ///< cell2D edge directions
        mesh_geometric_data_config.Cell2DsEdgesCentroid = true;  ///< cell2D edge centroid
        mesh_geometric_data_config.Cell2DsEdgeLengths = true;    ///< cell2D edge lengths
        mesh_geometric_data_config.Cell2DsEdgeTangents = true;   ///< cell2D edge tangents
        mesh_geometric_data_config.Cell2DsEdgeNormals = true;    ///< cell2D edge normals
        mesh_geometric_data_config.Cell2DsChebyshevCenter = false;
        mesh_geometric_data_config.Cell2DsTriangulationsByChebyshevCenter = false;            ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsInRadius = false; ///< cell2D triangulations
    }
};



struct FEM_PCC_2D_ReferenceElement final
{

    FEM_PCC_2D_ReferenceElement()
    {
    }
    ~FEM_PCC_2D_ReferenceElement(){};

    FEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {

        Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data result;

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;

        Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement triangle_reference_element;
        Polydim::FEM::PCC::FEM_Quadrilateral_PCC_2D_ReferenceElement quadrilateral_reference_element;

        result.triangle_reference_element_data = triangle_reference_element.Create(order);
        result.quadrilateral_reference_element_data = quadrilateral_reference_element.Create(order);

        return result;
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
