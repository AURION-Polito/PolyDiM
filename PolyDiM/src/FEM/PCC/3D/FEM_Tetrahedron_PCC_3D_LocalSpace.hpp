#ifndef __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP
#define __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP

#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"
#include "MapTetrahedron.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;
    double Tolerance3D;

    Eigen::MatrixXd Vertices;
};

struct FEM_Tetrahedron_PCC_3D_LocalSpace_Data final
{
    Gedim::MapTetrahedron::MapTetrahedronData MapData;
    unsigned int Order;                                   ///< Order of the space
    unsigned int NumberOfBasisFunctions;                  ///< Number of basis functions
    Eigen::MatrixXd Dofs;                                 ///< DOFs geometric position
    Gedim::Quadrature::QuadratureData InternalQuadrature; ///< Internal quadrature points and weights
};

/// \brief Interface used to FEM Values computation
class FEM_Tetrahedron_PCC_3D_LocalSpace final
{
  private:
    /// \brief map basis function values on element with correct order
    /// \note local_space not used for order 1
    inline Eigen::MatrixXd MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &, const Eigen::MatrixXd &referenceValues) const
    {
        return referenceValues;
    }

    /// \brief map basis function derivative values on element with correct order
    std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const;

  public:
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data CreateLocalSpace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                            const FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = Gedim::MapTetrahedron::FInv(local_space.MapData, points);

        FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                              const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = Gedim::MapTetrahedron::FInv(local_space.MapData, points);

        FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
