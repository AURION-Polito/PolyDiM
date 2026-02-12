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

#ifndef __I_ZFEM_PCC_2D_ReferenceElement_HPP
#define __I_ZFEM_PCC_2D_ReferenceElement_HPP

#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "MeshUtilities.hpp"
#include "Monomials_Data.hpp"

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{
struct ZFEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data fem_reference_element_data;
    Utilities::Monomials_Data monomials_data;

    Gedim::MeshUtilities::MeshGeometricData2DConfig mesh_geometric_data_config;

    ZFEM_PCC_2D_ReferenceElement_Data()
    {
        mesh_geometric_data_config.Cell2DsBoundingBox = false;
        mesh_geometric_data_config.Cell2DsTriangulations = false; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsAreas = true;          ///< cell2D areas
        mesh_geometric_data_config.Cell2DsCentroids = false;      ///< cell2D centroids
        mesh_geometric_data_config.Cell2DsDiameters = true;      ///< cell2D diameters
        mesh_geometric_data_config.Cell2DsEdgeDirections = true; ///< cell2D edge directions
        mesh_geometric_data_config.Cell2DsEdgesCentroid = false;  ///< cell2D edge centroid
        mesh_geometric_data_config.Cell2DsEdgeLengths = true;    ///< cell2D edge lengths
        mesh_geometric_data_config.Cell2DsEdgeTangents = true;   ///< cell2D edge tangents
        mesh_geometric_data_config.Cell2DsEdgeNormals = true;    ///< cell2D edge normals
        mesh_geometric_data_config.Cell2DsChebyshevCenter = true;
        mesh_geometric_data_config.Cell2DsTriangulationsByChebyshevCenter = true;            ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsInRadius = true; ///< cell2D triangulations
    }
};

class I_ZFEM_PCC_2D_ReferenceElement
{
  public:
    virtual ZFEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
