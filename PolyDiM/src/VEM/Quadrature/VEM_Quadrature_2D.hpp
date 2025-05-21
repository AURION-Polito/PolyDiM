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

#ifndef __VEM_Quadrature_2D_HPP
#define __VEM_Quadrature_2D_HPP

#include "Eigen/Eigen"
#include "QuadratureData.hpp"

namespace Polydim
{
namespace VEM
{
namespace Quadrature
{
struct VEM_QuadratureData_2D
{
    Gedim::Quadrature::QuadratureData ReferenceSegmentQuadrature;
    Gedim::Quadrature::QuadratureData ReferenceEdgeDOFsQuadrature;
    Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;

    Eigen::MatrixXd ReferenceSegmentInternalPoints;
    Eigen::VectorXd ReferenceSegmentInternalWeights;
    Eigen::Vector2d ReferenceSegmentExtremaWeights;

    Eigen::MatrixXd ReferenceEdgeDOFsInternalPoints;
    Eigen::VectorXd ReferenceEdgeDOFsInternalWeights;
    Eigen::Vector2d ReferenceEdgeDOFsExtremaWeights;
};

class VEM_Quadrature_2D final
{
  public:
    struct Edges_QuadratureData
    {
        Gedim::Quadrature::QuadratureData Quadrature;
        std::vector<Eigen::VectorXd> WeightsTimesNormal;
    };

    VEM_QuadratureData_2D Compute_PCC_2D(const unsigned int order) const;
    VEM_QuadratureData_2D Compute_MCC_2D(const unsigned int order) const;
    VEM_QuadratureData_2D Compute_MCC_EdgeOrtho_2D(const unsigned int order) const;
    VEM_QuadratureData_2D Compute_DF_PCC_2D(const unsigned int order) const;
    VEM_QuadratureData_2D Compute_DF_PCC_3D(const unsigned int order) const;

    Gedim::Quadrature::QuadratureData PolygonInternalQuadrature(const Gedim::Quadrature::QuadratureData &data,
                                                                const std::vector<Eigen::Matrix3d> &polygonTriangulationVertices) const;

    Edges_QuadratureData PolygonEdgesLobattoQuadrature(const Eigen::MatrixXd &ReferenceSegmentInternalPoints,
                                                       const Eigen::VectorXd &ReferenceSegmentInternalWeights,
                                                       const Eigen::Vector2d &ReferenceSegmentExtremaWeights,
                                                       const Eigen::MatrixXd &polygonVertices,
                                                       const Eigen::VectorXd &edgeLengths,
                                                       const std::vector<bool> &edgeDirections,
                                                       const Eigen::MatrixXd &edgeTangents,
                                                       const Eigen::MatrixXd &edgeNormals) const;

    VEM_Quadrature_2D::Edges_QuadratureData PolygonEdgesQuadrature(const Gedim::Quadrature::QuadratureData &data,
                                                                   const Eigen::MatrixXd &polygonVertices,
                                                                   const Eigen::VectorXd &edgeLengths,
                                                                   const std::vector<bool> &edgeDirections,
                                                                   const Eigen::MatrixXd &edgeTangents,
                                                                   const Eigen::MatrixXd &edgeNormals) const;
};
} // namespace Quadrature
} // namespace VEM
} // namespace Polydim

#endif
