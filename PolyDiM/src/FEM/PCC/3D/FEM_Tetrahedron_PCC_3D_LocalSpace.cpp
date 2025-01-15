#include "FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_Tetrahedron_PCC_3D_LocalSpace_Data FEM_Tetrahedron_PCC_3D_LocalSpace::CreateLocalSpace(
    const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
    const FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = polyhedron.Tolerance1D;
    geometry_utilities_config.Tolerance2D = polyhedron.Tolerance2D;
    geometry_utilities_config.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MapTetrahedron mapTetrehedron(geometry_utilities);
    localSpace.MapData = mapTetrehedron.Compute(polyhedron.Vertices);

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    for (unsigned int e = 0; e < 6; ++e)
    {
        const unsigned int edge_origin_index = polyhedron.Edges(0, e);
        const unsigned int edge_end_index = polyhedron.Edges(1, e);

        localSpace.polyhedron_to_reference_edge_index[e] =
            reference_element_data.Edges_by_vertices.at({edge_origin_index, edge_end_index});
    }

    for (unsigned int f = 0; f < 4; ++f)
    {
        const unsigned int face_edge_index = polyhedron.Faces[f](1, 0);
        const unsigned int face_reference_edge_index = localSpace.polyhedron_to_reference_edge_index[face_edge_index];
        const unsigned int face_vertex_index = polyhedron.Faces[f](0, 2);

        localSpace.polyhedron_to_reference_face_index[f] =
            reference_element_data.Faces_by_edge_vertex.at({face_reference_edge_index, face_vertex_index});
    }

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);

    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 4; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        unsigned int vertex_dofCounter = reference_element_data.NumDofs0D * v;
        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[vertex_dofCounter] = d;
            vertex_dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[4]);
    for (unsigned int e = 0; e < 6; ++e)
    {
        localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

        const unsigned int ref_e = localSpace.polyhedron_to_reference_edge_index[e];
        unsigned int edge_dof_counter = reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * ref_e;

        if (polyhedron.EdgesDirection.at(e))
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e]; d < localSpace.Dof1DsIndex[e + 1]; d++)
            {
                localSpace.DofsMeshOrder[edge_dof_counter] = d;
                edge_dof_counter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e + 1] - 1; d < UINT_MAX && d >= localSpace.Dof1DsIndex[e]; d--)
            {
                localSpace.DofsMeshOrder[edge_dof_counter] = d;
                edge_dof_counter++;
            }
        }
    }

    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[6]);
    for (unsigned int f = 0; f < 4; ++f)
    {
        localSpace.Dof2DsIndex[f + 1] = localSpace.Dof2DsIndex[f] + reference_element_data.NumDofs2D;

        const unsigned int ref_f = localSpace.polyhedron_to_reference_face_index[f];
        unsigned int face_dof_counter = reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * 6 +
                                        reference_element_data.NumDofs2D * ref_f;
        if (polyhedron.FacesDirection.at(f))
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
            {
                localSpace.DofsMeshOrder[face_dof_counter] = d;
                face_dof_counter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f]; d--)
            {
                localSpace.DofsMeshOrder[face_dof_counter] = d;
                face_dof_counter++;
            }
        }
    }

    localSpace.Dof3DsIndex.fill(localSpace.Dof2DsIndex[4]);
    localSpace.Dof3DsIndex[1] = localSpace.Dof3DsIndex[0] + reference_element_data.NumDofs3D;

    unsigned int cell_dof_counter =
        reference_element_data.NumDofs0D * 4 + reference_element_data.NumDofs1D * 6 + reference_element_data.NumDofs2D * 4;
    for (unsigned int d = localSpace.Dof3DsIndex[0]; d < localSpace.Dof3DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[cell_dof_counter] = d;
        cell_dof_counter++;
    }

    localSpace.Dofs = MapValues(localSpace, Gedim::MapTetrahedron::F(localSpace.MapData, reference_element_data.DofPositions));

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceTetrahedronQuadrature, localSpace.MapData);
    localSpace.BoundaryQuadrature =
        BoundaryQuadrature(reference_element_data.BoundaryReferenceElement_Data.ReferenceTriangleQuadrature, polyhedron);

    return localSpace;
}
// ***************************************************************************
MatrixXd FEM_Tetrahedron_PCC_3D_LocalSpace::MapValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                      const Eigen::MatrixXd &referenceValues) const
{
    Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
        basisFunctionValuesOrdered.col(local_space.DofsMeshOrder.at(d)) << referenceValues.col(d);

    return basisFunctionValuesOrdered;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_Tetrahedron_PCC_3D_LocalSpace::MapDerivativeValues(const FEM_Tetrahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                             const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
{
    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
        3,
        Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int i = 0; i < 3; i++)
    {
        basisFunctionsDerivativeValues[i] = local_space.MapData.QInv(i, i) * referenceDerivateValues[i];
        for (unsigned int j = 0; j < i; j++)
        {
            basisFunctionsDerivativeValues[i] += local_space.MapData.QInv(j, i) * referenceDerivateValues[j];
            basisFunctionsDerivativeValues[j] += local_space.MapData.QInv(i, j) * referenceDerivateValues[i];
        }
    }

    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
        3,
        Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
    {
        const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
        basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
        basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
        basisFunctionsDerivativeValuesOrdered.at(2).col(dofOrder) << basisFunctionsDerivativeValues.at(2).col(d);
    }

    return basisFunctionsDerivativeValuesOrdered;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Tetrahedron_PCC_3D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const Gedim::MapTetrahedron::MapTetrahedronData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    quadrature.Points = Gedim::MapTetrahedron::F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapTetrahedron::DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
std::vector<Gedim::Quadrature::QuadratureData> FEM_Tetrahedron_PCC_3D_LocalSpace::BoundaryQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const FEM_Tetrahedron_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    const unsigned int num_faces = polyhedron.Faces.size();
    std::vector<Gedim::Quadrature::QuadratureData> faces_quadrature(num_faces);

    for (unsigned int f = 0; f < num_faces; ++f)
    {
        auto &face_quadrature = faces_quadrature.at(f);

        const double face_area = polyhedron.FacesArea[f];
        const Eigen::Vector3d &face_translation = polyhedron.FacesTranslation[f];
        const Eigen::Matrix3d &face_rotation = polyhedron.FacesRotationMatrix[f];

        face_quadrature.Points = (face_rotation * reference_quadrature.Points).colwise() + face_translation;
        face_quadrature.Weights = reference_quadrature.Weights * std::abs(face_area);
    }

    return faces_quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
