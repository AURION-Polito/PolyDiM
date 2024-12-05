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
          VEM_QuadratureData_2D Quadrature2DData;
          Gedim::Quadrature::QuadratureData ReferenceSegmentQuadrature;
          Gedim::Quadrature::QuadratureData ReferenceTetrahedronQuadrature;
          Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;
          Eigen::MatrixXd ReferenceSegmentInternalPoints;
          Eigen::VectorXd ReferenceSegmentInternalWeights;
          Eigen::Vector2d ReferenceSegmentExtremaWeights;
      };

      class VEM_Quadrature_3D final
      {
        private:
          VEM_Quadrature_2D quadrature2D;

        public:
          struct Faces_QuadratureData
          {
              Gedim::Quadrature::QuadratureData Quadrature;
              std::vector<Eigen::VectorXd> WeightsTimesNormal;
          };

          struct Faces_QuadratureData2D
          {
              std::vector<Gedim::Quadrature::QuadratureData> FacesQuadrature;
              Gedim::Quadrature::QuadratureData Quadrature;
          };

          VEM_QuadratureData_3D Compute_PCC_3D(const unsigned int order) const;

          Gedim::Quadrature::QuadratureData PolyhedronInternalQuadrature(const VEM_QuadratureData_3D& data,
                                                                         const Gedim::GeometryUtilities& geometryUtility,
                                                                         const std::vector<Eigen::MatrixXd>& polyhedronTetrahedronVertices) const;

          Faces_QuadratureData PolyhedronFacesQuadrature(const Gedim::GeometryUtilities& geometryUtility,
                                                         const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                         const std::vector<Eigen::Matrix3d>& facesRotationMatrix,
                                                         const std::vector<Eigen::Vector3d>& facesTranslation,
                                                         const std::vector<Eigen::Vector3d>& facesNormals,
                                                         const std::vector<bool>& faceNormalDirections,
                                                         const std::vector<Eigen::MatrixXd>& facesQuadraturePoints,
                                                         const std::vector<Eigen::VectorXd>& facesQuadratureWeights) const;

          Faces_QuadratureData2D PolyhedronFacesQuadrature(const VEM_QuadratureData_3D& data,
                                                           const Gedim::GeometryUtilities& geometryUtility,
                                                           const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                                           const std::vector<std::vector<Eigen::Matrix3d> >& facesTriangulations2D,
                                                           const std::vector<Eigen::Matrix3d>& facesRotationMatrix,
                                                           const std::vector<Eigen::Vector3d>& facesTranslation) const;

          Eigen::MatrixXd PolyhedronInternalEdgesQuadraturePoints(const VEM_QuadratureData_3D& data,
                                                                  const Eigen::MatrixXd& polyhedronVertices,
                                                                  const Eigen::MatrixXi& polyhedronEdges,
                                                                  const std::vector<bool>& edgeDirections,
                                                                  const Eigen::MatrixXd& edgeTangents) const;

      };
    }
  }
}

#endif
