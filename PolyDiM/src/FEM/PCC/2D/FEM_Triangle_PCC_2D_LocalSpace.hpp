#ifndef __FEM_Triangle_PCC_2D_LocalSpace_HPP
#define __FEM_Triangle_PCC_2D_LocalSpace_HPP

#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "MapTriangle.hpp"

namespace Polydim
{
  namespace FEM
  {
    namespace PCC
    {
      struct FEM_Triangle_PCC_2D_Polygon_Geometry final
      {
          double Tolerance1D;
          double Tolerance2D;

          Eigen::MatrixXd Vertices;
          std::vector<bool> EdgesDirection;
          Eigen::MatrixXd EdgesTangent;
          Eigen::VectorXd EdgesLength;
      };

      struct FEM_Triangle_PCC_2D_LocalSpace_Data final
      {
          Gedim::MapTriangle::MapTriangleData MapData;
          Eigen::Matrix3d B_lap;
          unsigned int Order;                                   ///< Order of the space
          unsigned int NumberOfBasisFunctions;                  ///< Number of basis functions
          Eigen::MatrixXd Dofs;                                 ///< DOFs geometric position
          std::vector<unsigned int> DofsMeshOrder;              ///< DOFs position depending on element
          std::array<unsigned int, 4> Dof0DsIndex;              ///< local DOF index for each element 0D
          std::array<unsigned int, 4> Dof1DsIndex;              ///< local DOF index for each element 1D
          std::array<unsigned int, 2> Dof2DsIndex;              ///< local DOF index for each element 2D
          Gedim::Quadrature::QuadratureData InternalQuadrature; ///< Internal quadrature points and weights
          std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature; ///< Boundary quadrature points and weights on each edge
      };

      /// \brief Interface used to FEM Values computation
      class FEM_Triangle_PCC_2D_LocalSpace final
      {
        private:
          /// \brief map basis function values on element with correct order
          Eigen::MatrixXd MapValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

          /// \brief map basis function derivative values on element with correct order
          std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                           const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

          /// \brief map basis function derivative values on element with correct order
          Eigen::MatrixXd MapLaplacianValues(const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                             const std::array<Eigen::MatrixXd, 4> &referenceSecondDerivateValues) const;

          Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                               const Gedim::MapTriangle::MapTriangleData &mapData) const;

          std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                            const FEM_Triangle_PCC_2D_Polygon_Geometry &polygon) const;

        public:
          FEM_Triangle_PCC_2D_LocalSpace_Data CreateLocalSpace(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                               const FEM_Triangle_PCC_2D_Polygon_Geometry &polygon) const;

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
          {
            return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                      const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
          {
            if (local_space.Order > 2)
              throw std::runtime_error("Unsupported order");

            return MapLaplacianValues(local_space, reference_element_data.ReferenceBasisFunctionSecondDerivativeValues);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                                    const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space) const
          {
            return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                             const Eigen::MatrixXd &points) const
          {
            Gedim::MapTriangle mapTriangle;
            const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

            FEM_Triangle_PCC_2D_ReferenceElement reference_element;

            return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                                    const FEM_Triangle_PCC_2D_LocalSpace_Data &local_space,
                                                                                    const Eigen::MatrixXd &points) const
          {
            Gedim::MapTriangle mapTriangle;
            const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

            FEM_Triangle_PCC_2D_ReferenceElement reference_element;

            return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnEdge(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data) const
          {
            return reference_element_data.BoundaryReferenceElement_Data.ReferenceBasisFunctionValues;
          }
      };
    } // namespace PCC
  } // namespace FEM
} // namespace Polydim

#endif
