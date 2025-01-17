#include "VEM_DF_PCC_3D_Pressure_LocalSpace.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
//****************************************************************************
VEM_DF_PCC_3D_Pressure_LocalSpace_Data VEM_DF_PCC_3D_Pressure_LocalSpace::CreateLocalSpace(
    const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data_3D,
    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    VEM_DF_PCC_3D_Pressure_LocalSpace_Data localSpace;
    Quadrature::VEM_Quadrature_3D quadrature3D;

    localSpace.InternalQuadrature =
        quadrature3D.PolyhedronInternalQuadrature(reference_element_data_3D.Quadrature, geometryUtilities, polyhedron.TetrahedronVertices);

    InitializeProjectorsComputation(reference_element_data_3D,
                                    polyhedron.Centroid,
                                    polyhedron.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_DF_PCC_3D_Pressure_LocalSpace::InitializeProjectorsComputation(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                        const Eigen::Vector3d &polyhedronCentroid,
                                                                        const double &polyhedronDiameter,
                                                                        const Eigen::MatrixXd &internalQuadraturePoints,
                                                                        const Eigen::VectorXd &internalQuadratureWeights,
                                                                        VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const
{

    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumBasisFunctions = reference_element_data.NumDofs3D;

    localSpace.Nk =
        (reference_element_data.Order + 1) * (reference_element_data.Order + 2) * (reference_element_data.Order + 3) / 6;
    localSpace.Nkm1 = reference_element_data.Order * (reference_element_data.Order + 1) * (reference_element_data.Order + 2) / 6;

    // Compute Vandermonde matrices.
    localSpace.Diameter = polyhedronDiameter;
    localSpace.Centroid = polyhedronCentroid;

    localSpace.VanderInternal =
        monomials.Vander(reference_element_data.Monomials, internalQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;
}
//****************************************************************************
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim
