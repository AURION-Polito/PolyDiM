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

#ifndef __VEM_Quadrature_3D_HPP
#define __VEM_Quadrature_3D_HPP

#include "GeometryUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace Quadrature
{
struct VEM_QuadratureData_3D final
{
    Gedim::Quadrature::QuadratureData ReferenceTetrahedronQuadrature;
    Polydim::VEM::Quadrature::VEM_QuadratureData_2D QuadratureData_2D;
};

class VEM_Quadrature_3D final
{
  public:
    struct Faces_QuadratureData_PCC
    {
        Gedim::Quadrature::QuadratureData Quadrature;
        std::vector<Eigen::VectorXd> WeightsTimesNormal;
    };

    struct Faces_QuadratureData_MCC
    {
        Gedim::Quadrature::QuadratureData Quadrature;
        std::vector<Gedim::Quadrature::QuadratureData> FacesQuadrature;
    };

    Polydim::VEM::Quadrature::VEM_QuadratureData_3D Compute_PCC_3D(const unsigned int order) const;
    Polydim::VEM::Quadrature::VEM_QuadratureData_3D Compute_MCC_3D(const unsigned int order) const;
    Polydim::VEM::Quadrature::VEM_QuadratureData_3D Compute_DF_PCC_3D(const unsigned int order) const;

    Gedim::Quadrature::QuadratureData PolyhedronInternalQuadrature(const Polydim::VEM::Quadrature::VEM_QuadratureData_3D &data,
                                                                   const Gedim::GeometryUtilities &geometryUtility,
                                                                   const std::vector<Eigen::MatrixXd> &polyhedronTetrahedronVertices) const;

    Gedim::Quadrature::QuadratureData PolyhedronInternalQuadrature(const Gedim::GeometryUtilities &geometryUtility,
                                                                   const Gedim::Quadrature::QuadratureData &data,
                                                                   const std::vector<Eigen::MatrixXd> &polyhedronTetrahedronVertices) const;

    Polydim::VEM::Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC PolyhedronFacesQuadrature(
        const Gedim::GeometryUtilities &geometryUtility,
        const std::vector<Eigen::MatrixXi> &polyhedronFaces,
        const std::vector<Eigen::Matrix3d> &facesRotationMatrix,
        const std::vector<Eigen::Vector3d> &facesTranslation,
        const std::vector<Eigen::Vector3d> &facesNormals,
        const std::vector<bool> &faceNormalDirections,
        const std::vector<Eigen::MatrixXd> &facesQuadraturePoints,
        const std::vector<Eigen::VectorXd> &facesQuadratureWeights) const;

    Eigen::MatrixXd PolyhedronInternalEdgesQuadraturePoints(const Eigen::MatrixXd &referenceSegmentInternalPoints,
                                                            const Eigen::MatrixXd &polyhedronVertices,
                                                            const Eigen::MatrixXi &polyhedronEdges,
                                                            const std::vector<bool> &edgeDirections,
                                                            const Eigen::MatrixXd &edgeTangents) const;

    Polydim::VEM::Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_MCC PolyhedronFacesQuadrature(
        const Polydim::VEM::Quadrature::VEM_QuadratureData_3D &data,
        const Gedim::GeometryUtilities &geometryUtility,
        const std::vector<std::vector<Eigen::Matrix3d>> &facesTriangulations2D,
        const std::vector<Eigen::Matrix3d> &facesRotationMatrix,
        const std::vector<Eigen::Vector3d> &facesTranslation) const;
};
} // namespace Quadrature
} // namespace VEM
} // namespace Polydim

#endif
