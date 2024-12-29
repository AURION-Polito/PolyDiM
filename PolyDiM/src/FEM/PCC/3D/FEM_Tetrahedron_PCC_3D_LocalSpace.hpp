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
    const double Tolerance1D;
    const double Tolerance2D;
    const double Tolerance3D;

    const Eigen::MatrixXd Vertices;
};

struct FEM_Tetrahedron_PCC_3D_LocalSpace_Data final
{
    const double Tolerance1D;
    const double Tolerance2D;
    const double Tolerance3D;
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
    Eigen::MatrixXd MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

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
      Gedim::GeometryUtilitiesConfig geometry_utilities_config;
      geometry_utilities_config.Tolerance1D = local_space.Tolerance1D;
      geometry_utilities_config.Tolerance2D = local_space.Tolerance2D;
      geometry_utilities_config.Tolerance3D = local_space.Tolerance3D;
      Gedim::GeometryUtilities geometry_utilties(geometry_utilities_config);

        Gedim::MapTetrahedron mapTetrahedron(geometry_utilties);
        const Eigen::MatrixXd referencePoints = mapTetrahedron.FInv(local_space.MapData, points);

        FEM_RefElement_Langrange_PCC_Tetrahedron_3D reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                              const Eigen::MatrixXd &points) const
    {
      Gedim::GeometryUtilitiesConfig geometry_utilities_config;
      geometry_utilities_config.Tolerance1D = local_space.Tolerance1D;
      geometry_utilities_config.Tolerance2D = local_space.Tolerance2D;
      geometry_utilities_config.Tolerance3D = local_space.Tolerance3D;
      Gedim::GeometryUtilities geometry_utilties(geometry_utilities_config);

        Gedim::MapTetrahedron mapTetrahedron(geometry_utilties);
        const Eigen::MatrixXd referencePoints = mapTetrahedron.FInv(local_space.MapData, points);

        FEM_RefElement_Langrange_PCC_Tetrahedron_3D reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
