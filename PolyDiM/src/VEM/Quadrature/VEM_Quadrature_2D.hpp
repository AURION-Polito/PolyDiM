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
    Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;
    Eigen::MatrixXd ReferenceSegmentInternalPoints;
    Eigen::VectorXd ReferenceSegmentInternalWeights;
    Eigen::Vector2d ReferenceSegmentExtremaWeights;
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

    Gedim::Quadrature::QuadratureData PolygonInternalQuadrature(const VEM_QuadratureData_2D& data,
                                                                const std::vector<Eigen::Matrix3d>& polygonTriangulationVertices) const;

    Edges_QuadratureData PolygonEdgesLobattoQuadrature(const VEM_QuadratureData_2D& data,
                                                       const Eigen::MatrixXd& polygonVertices,
                                                       const Eigen::VectorXd& edgeLengths,
                                                       const std::vector<bool>& edgeDirections,
                                                       const Eigen::MatrixXd& edgeTangents,
                                                       const Eigen::MatrixXd& edgeNormals) const;

    VEM_Quadrature_2D::Edges_QuadratureData PolygonEdgesQuadrature(const VEM_QuadratureData_2D &data,
                                                                   const Eigen::MatrixXd &polygonVertices,
                                                                   const Eigen::VectorXd &edgeLengths,
                                                                   const std::vector<bool> &edgeDirections,
                                                                   const Eigen::MatrixXd &edgeTangents,
                                                                   const Eigen::MatrixXd &edgeNormals) const;

};
}
}
}

#endif
