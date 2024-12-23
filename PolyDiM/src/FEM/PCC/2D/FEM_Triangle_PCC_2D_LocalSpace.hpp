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
          const double Tolerance1D;
          const double Tolerance2D;

          const Eigen::MatrixXd &Vertices;
          const std::vector<bool> &EdgesDirection;
      };

      struct FEM_Triangle_PCC_2D_LocalSpace_Data final
      {
          Gedim::MapTriangle::MapTriangleData MapData;
          unsigned int Order; ///< Order of the space
          unsigned int NumberOfBasisFunctions; ///< Number of basis functions
          Eigen::MatrixXd Dofs; ///< DOFs geometric position
          std::vector<unsigned int> DofsMeshOrder; ///< DOFs position depending on mesh directions
          std::vector<unsigned int> Dof0DsIndex; ///< local DOF index for each element 0D, size num0D + 1
          std::vector<unsigned int> Dof1DsIndex; ///< local DOF index for each element 1D, size num1D + 1
          std::vector<unsigned int> Dof2DsIndex; ///< local DOF index for each element 2D, size num2D + 1
          std::vector<unsigned int> Dof3DsIndex; ///< local DOF index for each element 3D, size num3D + 1
          Gedim::Quadrature::QuadratureData InternalQuadrature; ///< Internal quadrature points and weights
      };

      /// \brief Interface used to FEM Values computation
      class FEM_Triangle_PCC_2D_LocalSpace final
      {
        private:
          /// \brief map basis function values on element with correct order
          Eigen::MatrixXd MapValues(const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space,
                                    const Eigen::MatrixXd& referenceValues) const;

          /// \brief map basis function derivative values on element with correct order
          std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space,
                                                           const std::vector<Eigen::MatrixXd>& referenceDerivateValues) const;

          Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData& reference_quadrature,
                                                               const Gedim::MapTriangle::MapTriangleData& mapData) const;

        public:
          FEM_Triangle_PCC_2D_LocalSpace_Data CreateLocalSpace(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                               const FEM_Triangle_PCC_2D_Polygon_Geometry &polygon) const;

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space) const
          {
            return MapValues(local_space,
                             reference_element_data.ReferenceBasisFunctionValues);
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                                    const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space) const
          {
            return MapDerivativeValues(local_space,
                                       reference_element_data.ReferenceBasisFunctionDerivativeValues);
          }

          inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space,
                                                             const Eigen::MatrixXd& points) const
          {
            Gedim::MapTriangle mapTriangle;
            const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData,
                                                                     points);

            FEM_RefElement_Langrange_PCC_Triangle_2D reference_element;

            return MapValues(local_space,
                             reference_element.EvaluateBasisFunctions(referencePoints,
                                                                      reference_element_data));
          }

          inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                                    const FEM_Triangle_PCC_2D_LocalSpace_Data& local_space,
                                                                                    const Eigen::MatrixXd& points) const
          {
            Gedim::MapTriangle mapTriangle;
            const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData,
                                                                     points);

            FEM_RefElement_Langrange_PCC_Triangle_2D reference_element;

            return MapDerivativeValues(local_space,
                                       reference_element.EvaluateBasisFunctionDerivatives(referencePoints,
                                                                                          reference_element_data));
          }
      };
    }
  }
}

#endif
