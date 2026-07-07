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

/// @brief Reference-element data for the 2D primal conforming ZFEM (Zipped FEM) space.
///
/// Holds the order- and dimension-dependent quantities defining the ZFEM reference
/// element: the DOF counts per mesh-entity dimension, the underlying triangular FEM
/// reference element, the monomial basis data, and the mesh geometric-data
/// configuration required by the method. The default constructor enables exactly the
/// geometric quantities ZFEM needs on each polygon (areas, diameters, edge
/// directions/lengths/tangents/normals, Chebyshev center and its triangulation, and
/// in-radius), disabling the others to avoid unnecessary computation.
struct ZFEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension; ///< Spatial dimension of the element (2).
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of DOFs per vertex (0D entity).
    unsigned int NumDofs1D; ///< Number of DOFs per edge (1D entity).
    unsigned int NumDofs2D; ///< Number of DOFs internal to the cell (2D entity).

    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data fem_reference_element_data; ///< Underlying triangular
                                                                                             ///< FEM reference element.
    Utilities::Monomials_Data monomials_data;                                                ///< Monomial basis data.

    Gedim::MeshUtilities::MeshGeometricData2DConfig mesh_geometric_data_config; ///< Geometric quantities required by
                                                                                ///< the method.

    /// @brief Default constructor: configure the required mesh geometric quantities.
    ///
    /// Enables the geometric data ZFEM depends on and disables the rest.
    ZFEM_PCC_2D_ReferenceElement_Data()
    {
        mesh_geometric_data_config.Cell2DsBoundingBox = false;
        mesh_geometric_data_config.Cell2DsTriangulations = false; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsAreas = true;           ///< cell2D areas
        mesh_geometric_data_config.Cell2DsCentroids = false;      ///< cell2D centroids
        mesh_geometric_data_config.Cell2DsDiameters = true;       ///< cell2D diameters
        mesh_geometric_data_config.Cell2DsEdgeDirections = true;  ///< cell2D edge directions
        mesh_geometric_data_config.Cell2DsEdgesCentroid = false;  ///< cell2D edge centroid
        mesh_geometric_data_config.Cell2DsEdgeLengths = true;     ///< cell2D edge lengths
        mesh_geometric_data_config.Cell2DsEdgeTangents = true;    ///< cell2D edge tangents
        mesh_geometric_data_config.Cell2DsEdgeNormals = true;     ///< cell2D edge normals
        mesh_geometric_data_config.Cell2DsChebyshevCenter = true;
        mesh_geometric_data_config.Cell2DsTriangulationsByChebyshevCenter = true; ///< cell2D triangulations
        mesh_geometric_data_config.Cell2DsInRadius = true;                        ///< cell2D triangulations
    }
};

/// @brief Abstract interface for the 2D primal conforming ZFEM reference element.
///
/// Defines the factory contract that concrete ZFEM reference-element implementations
/// must fulfill, producing a fully initialized @ref ZFEM_PCC_2D_ReferenceElement_Data
/// for a requested polynomial order.
class I_ZFEM_PCC_2D_ReferenceElement
{
  public:
    virtual ZFEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
