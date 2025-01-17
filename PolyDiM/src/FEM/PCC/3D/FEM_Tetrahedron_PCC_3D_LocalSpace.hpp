#ifndef __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP
#define __FEM_Tetrahedron_PCC_3D_LocalSpace_HPP

#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"
#include "MapTetrahedron.hpp"
#include "MapTriangle.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_Tetrahedron_PCC_3D_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;
    double Tolerance3D;

    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    std::vector<Eigen::MatrixXd> Faces_2D_Vertices;
    std::vector<bool> EdgesDirection;
    std::vector<bool> FacesDirection;
    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
};

struct FEM_Tetrahedron_PCC_3D_LocalSpace_Data final
{
    Gedim::MapTetrahedron::MapTetrahedronData MapData;
    std::array<Gedim::MapTriangle::MapTriangleData, 4> Faces_MapData;
    unsigned int Order;                  ///< Order of the space
    unsigned int NumberOfBasisFunctions; ///< Number of basis functions
    Eigen::MatrixXd Dofs;                ///< DOFs geometric position
    std::array<unsigned int, 6> polyhedron_to_reference_edge_index;
    std::array<bool, 6> polyhedron_to_reference_edge_direction;
    std::array<unsigned int, 4> polyhedron_to_reference_face_index;
    std::vector<unsigned int> DofsMeshOrder;                           ///< DOFs position depending on element
    std::array<unsigned int, 5> Dof0DsIndex;                           ///< local DOF index for each element 0D
    std::array<unsigned int, 7> Dof1DsIndex;                           ///< local DOF index for each element 1D
    std::array<unsigned int, 5> Dof2DsIndex;                           ///< local DOF index for each element 2D
    std::array<unsigned int, 2> Dof3DsIndex;                           ///< local DOF index for each element 3D
    Gedim::Quadrature::QuadratureData InternalQuadrature;              ///< Internal quadrature points and weights
    std::array<Gedim::Quadrature::QuadratureData, 4> BoundaryQuadrature; ///< Boundary quadrature points and weights on
                                                                       ///< each face
};

/// \brief Interface used to FEM Values computation
class FEM_Tetrahedron_PCC_3D_LocalSpace final
{
  private:
    /// \brief map basis function values on element with correct order
    /// \note local_space not used for order 1
    Eigen::MatrixXd MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space, const Eigen::MatrixXd &referenceValues) const;

    /// \brief map basis function derivative values on element with correct order
    std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const;
    std::array<Gedim::Quadrature::QuadratureData, 4> BoundaryQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
        const std::array<Gedim::MapTriangle::MapTriangleData, 4> &faces_mapData,
        const FEM_Tetrahedron_PCC_3D_Geometry &polyhedron) const;

  public:
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data CreateLocalSpace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                            const FEM_Tetrahedron_PCC_3D_Geometry &polyhedron) const;

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

    inline Eigen::MatrixXd ComputeBasisFunctionsValuesOnFace(const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {
        return reference_element_data.BoundaryReferenceElement_Data.ReferenceBasisFunctionValues;
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
