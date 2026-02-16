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

#ifndef __I_VEM_PCC_2D_ReferenceElement_HPP
#define __I_VEM_PCC_2D_ReferenceElement_HPP

#include "MeshUtilities.hpp"
#include "Monomials_Data.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
struct VEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    Utilities::Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_2D Quadrature;

    Gedim::MeshUtilities::MeshGeometricData2DConfig mesh_geometric_data_config;

    VEM_PCC_2D_ReferenceElement_Data()
    {
        mesh_geometric_data_config.Cell2DsBoundingBox = false;
        mesh_geometric_data_config.Cell2DsTriangulations = true; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsAreas = true;          ///< cell2D areas
        mesh_geometric_data_config.Cell2DsCentroids = true;      ///< cell2D centroids
        mesh_geometric_data_config.Cell2DsDiameters = true;      ///< cell2D diameters
        mesh_geometric_data_config.Cell2DsEdgeDirections = true; ///< cell2D edge directions
        mesh_geometric_data_config.Cell2DsEdgesCentroid = false; ///< cell2D edge centroid
        mesh_geometric_data_config.Cell2DsEdgeLengths = true;    ///< cell2D edge lengths
        mesh_geometric_data_config.Cell2DsEdgeTangents = true;   ///< cell2D edge tangents
        mesh_geometric_data_config.Cell2DsEdgeNormals = true;    ///< cell2D edge normals
        mesh_geometric_data_config.Cell2DsChebyshevCenter = false;
        mesh_geometric_data_config.Cell2DsTriangulationsByChebyshevCenter = false; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsInRadius = false;                        ///< cell2D triangulations
    }
};

class I_VEM_PCC_2D_ReferenceElement
{
  public:
    virtual VEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
