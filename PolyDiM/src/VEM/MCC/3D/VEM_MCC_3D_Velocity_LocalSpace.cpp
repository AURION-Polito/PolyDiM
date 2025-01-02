#include "VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "LAPACK_utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace MCC
{
//****************************************************************************
VEM_MCC_3D_Velocity_LocalSpace_Data VEM_MCC_3D_Velocity_LocalSpace::CreateLocalSpace(const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                                     const VEM_MCC_3D_Polyhedron_Geometry &polyhedron) const
{
    VEM_MCC_3D_Velocity_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Quadrature::VEM_Quadrature_3D quadrature;
    localSpace.InternalQuadrature =
        quadrature.PolyhedronInternalQuadrature(reference_element_data.Quadrature, geometryUtilities, polyhedron.TetrahedronVertices);

    localSpace.BoundaryQuadrature = quadrature.PolyhedronFacesQuadrature(reference_element_data.Quadrature,
                                                                         geometryUtilities,
                                                                         polyhedron.FacesTriangulationVertices2D,
                                                                         polyhedron.FacesRotationMatrix,
                                                                         polyhedron.FacesTranslation);

    InitializeProjectorsComputation(reference_element_data,
                                    polyhedron.FacesMeasure.size(),
                                    polyhedron.Centroid,
                                    polyhedron.Measure,
                                    polyhedron.Diameter,
                                    localSpace.InternalQuadrature.Points,
                                    localSpace.InternalQuadrature.Weights,
                                    localSpace.BoundaryQuadrature.Quadrature.Points,
                                    localSpace);

    Eigen::MatrixXd W2;
    Eigen::MatrixXd B2Nabla;
    ComputeValuesOnBoundary(reference_element_data,
                            polyhedron.FacesMeasure.size(),
                            polyhedron.FacesNormal,
                            polyhedron.FacesNormalDirection,
                            polyhedron.FacesGlobalNormalDirection,
                            polyhedron.FacesCentroid2D,
                            polyhedron.FacesMeasure,
                            polyhedron.FacesDiameter,
                            localSpace.BoundaryQuadrature.FacesQuadrature,
                            localSpace.BoundaryQuadrature.Quadrature.Weights,
                            W2,
                            B2Nabla,
                            localSpace);

    ComputeDivergenceCoefficients(polyhedron.Measure, W2, localSpace);

    ComputeL2Projectors(polyhedron.Measure, localSpace.InternalQuadrature.Weights, B2Nabla, localSpace);

    ComputePolynomialBasisDofs(polyhedron.Measure, localSpace);

    return localSpace;
}
//****************************************************************************
void VEM_MCC_3D_Velocity_LocalSpace::InitializeProjectorsComputation(const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                     const unsigned int &numFaces,
                                                                     const Eigen::Vector3d &polyhedronCentroid,
                                                                     const double &polyhedronMeasure,
                                                                     const double &polyhedronDiameter,
                                                                     const Eigen::MatrixXd &internalQuadraturePoints,
                                                                     const Eigen::VectorXd &internalQuadratureWeights,
                                                                     const Eigen::MatrixXd &boundaryQuadraturePoints,
                                                                     VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    localSpace.Dimension = reference_element_data.Dimension;
    localSpace.Order = reference_element_data.Order;

    localSpace.NumBoundaryBasisFunctions = numFaces * reference_element_data.NumDofs2D;

    localSpace.Nk =
        (reference_element_data.Order + 1) * (reference_element_data.Order + 2) * (reference_element_data.Order + 3) / 6;

    localSpace.NkNabla = localSpace.Nk + (reference_element_data.Order + 2) * (reference_element_data.Order + 3) / 2 - 1;

    localSpace.NumNablaInternalBasisFunctions = localSpace.Nk - 1;
    localSpace.NumBigOPlusInternalBasisFunctions = 3 * localSpace.Nk - localSpace.NkNabla;
    localSpace.NumInternalBasisFunctions = localSpace.NumNablaInternalBasisFunctions + localSpace.NumBigOPlusInternalBasisFunctions;

    localSpace.NumBasisFunctions = localSpace.NumBoundaryBasisFunctions + localSpace.NumInternalBasisFunctions;

    // Compute Vandermonde matrices.
    localSpace.QmatrixKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                               reference_element_data.MonomialsKp1.NumMonomials);
    localSpace.QmatrixInvKp1 = MatrixXd::Identity(reference_element_data.MonomialsKp1.NumMonomials,
                                                  reference_element_data.MonomialsKp1.NumMonomials);

    localSpace.Centroid = polyhedronCentroid;
    localSpace.Diameter = polyhedronDiameter;
    localSpace.Measure = polyhedronMeasure;

    localSpace.VanderInternalKp1 =
        monomials3D.Vander(reference_element_data.MonomialsKp1, internalQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderInternal = localSpace.VanderInternalKp1.leftCols(localSpace.Nk);

    localSpace.VanderBoundaryKp1 =
        monomials3D.Vander(reference_element_data.MonomialsKp1, boundaryQuadraturePoints, polyhedronCentroid, polyhedronDiameter);

    localSpace.VanderBoundary = localSpace.VanderBoundaryKp1.leftCols(localSpace.Nk);

    // Compute mass matrix of monomials.
    const VectorXd sqrtInternalQuadratureWeights = internalQuadratureWeights.array().sqrt();
    const MatrixXd temp = sqrtInternalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    localSpace.Hmatrix = temp.transpose() * temp;

    localSpace.TkNabla.resize(localSpace.NkNabla, 3 * localSpace.Nk);

    const MatrixXd dx =
        (1.0 / polyhedronDiameter) *
        monomials3D.D_x(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);
    const MatrixXd dy =
        (1.0 / polyhedronDiameter) *
        monomials3D.D_y(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);
    const MatrixXd dz =
        (1.0 / polyhedronDiameter) *
        monomials3D.D_z(reference_element_data.MonomialsKp1).bottomLeftCorner(localSpace.NkNabla, localSpace.Nk);

    localSpace.TkNabla << dx, dy, dz;
    MatrixXd V;
    VectorXd S;
    LAPACK_utilities::svd(localSpace.TkNabla, V, S);

    MatrixXd VanderInternalEnlarged(localSpace.Dimension * localSpace.Nk,
                                    localSpace.Dimension * internalQuadratureWeights.size());

    const MatrixXd zeros = MatrixXd::Zero(localSpace.VanderInternal.cols(), localSpace.VanderInternal.rows());

    VanderInternalEnlarged << localSpace.VanderInternal.transpose(), zeros, zeros, zeros,
        localSpace.VanderInternal.transpose(), zeros, zeros, zeros, localSpace.VanderInternal.transpose();

    const MatrixXd GkNablaVanderInternal = localSpace.TkNabla * VanderInternalEnlarged;

    localSpace.TkBigOPlus = V.transpose().rightCols(localSpace.Dimension * localSpace.Nk - localSpace.NkNabla).transpose();

    const MatrixXd GkBigOPlusVanderInternal = localSpace.TkBigOPlus * VanderInternalEnlarged;

    localSpace.GkVanderInternal.resize(localSpace.Dimension * localSpace.Nk,
                                       localSpace.Dimension * internalQuadratureWeights.size());

    localSpace.GkVanderInternal << GkNablaVanderInternal, GkBigOPlusVanderInternal;

    VectorXd internalWeights3k(localSpace.Dimension * internalQuadratureWeights.size());
    internalWeights3k << sqrtInternalQuadratureWeights, sqrtInternalQuadratureWeights, sqrtInternalQuadratureWeights;

    const MatrixXd tempG = internalWeights3k.asDiagonal() * localSpace.GkVanderInternal.transpose();
    localSpace.Gmatrix = tempG.transpose() * tempG;
}
//****************************************************************************
void VEM_MCC_3D_Velocity_LocalSpace::ComputeDivergenceCoefficients(const double &polyhedronMeasure,
                                                                   const Eigen::MatrixXd &W2,
                                                                   VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    MatrixXd W1 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);

    if (localSpace.Order > 0)
    {
        W1.block(1, localSpace.NumBoundaryBasisFunctions, localSpace.NumNablaInternalBasisFunctions, localSpace.NumNablaInternalBasisFunctions) =
            -polyhedronMeasure * Eigen::MatrixXd::Identity(localSpace.NumNablaInternalBasisFunctions,
                                                           localSpace.NumNablaInternalBasisFunctions);
    }

    localSpace.Wmatrix = W1 + W2;
    localSpace.Vmatrix = localSpace.Hmatrix.llt().solve(localSpace.Wmatrix);
}
//****************************************************************************
void VEM_MCC_3D_Velocity_LocalSpace::ComputeL2Projectors(const double &polyhedronMeasure,
                                                         const Eigen::VectorXd &internalQuadratureWeights,
                                                         const Eigen::MatrixXd &B2Nabla,
                                                         VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    const MatrixXd HHashtagMatrix = localSpace.VanderInternalKp1.rightCols(localSpace.NkNabla).transpose() *
                                    internalQuadratureWeights.asDiagonal() * localSpace.VanderInternal;
    const MatrixXd B1Nabla = -HHashtagMatrix * localSpace.Vmatrix;

    const MatrixXd BNabla = B1Nabla + B2Nabla;

    MatrixXd BBigOPlus = MatrixXd::Zero(localSpace.NumBigOPlusInternalBasisFunctions, localSpace.NumBasisFunctions);
    BBigOPlus.rightCols(localSpace.NumBigOPlusInternalBasisFunctions) =
        polyhedronMeasure * Eigen::MatrixXd::Identity(localSpace.NumBigOPlusInternalBasisFunctions,
                                                      localSpace.NumBigOPlusInternalBasisFunctions);

    localSpace.Bmatrix.resize(localSpace.Dimension * localSpace.Nk, localSpace.NumBasisFunctions);
    localSpace.Bmatrix << BNabla, BBigOPlus;

    localSpace.Pi0k = localSpace.Gmatrix.llt().solve(localSpace.Bmatrix);
}
//****************************************************************************
void VEM_MCC_3D_Velocity_LocalSpace::ComputeValuesOnBoundary(const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                             const unsigned int &numFaces,
                                                             const std::vector<Vector3d> &facesNormals,
                                                             const std::vector<bool> &facesNormalDirections,
                                                             const std::vector<bool> &faceNormalGlobalDirections,
                                                             const std::vector<Vector3d> &facesCentroids,
                                                             const std::vector<double> &facesAreas,
                                                             const std::vector<double> &facesDiameters,
                                                             const std::vector<Gedim::Quadrature::QuadratureData> &facesQuadrature,
                                                             const Eigen::VectorXd &boundaryQuadratureWeights,
                                                             MatrixXd &W2,
                                                             MatrixXd &B2Nabla,
                                                             VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
{
    const unsigned int numBoundaryQuadraturePoints = boundaryQuadratureWeights.size();

    MatrixXd VanderBoundaryEnlarged(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * numBoundaryQuadraturePoints);

    const MatrixXd zeros = MatrixXd::Zero(localSpace.VanderBoundary.cols(), localSpace.VanderBoundary.rows());

    VanderBoundaryEnlarged << localSpace.VanderBoundary.transpose(), zeros, zeros, zeros,
        localSpace.VanderBoundary.transpose(), zeros, zeros, zeros, localSpace.VanderBoundary.transpose();

    const MatrixXd GkNablaVanderBoundary = localSpace.TkNabla * VanderBoundaryEnlarged;

    const MatrixXd GkBigOPlusVanderBoundary = localSpace.TkBigOPlus * VanderBoundaryEnlarged;

    MatrixXd GkVanderBoundary(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * numBoundaryQuadraturePoints);

    GkVanderBoundary << GkNablaVanderBoundary, GkBigOPlusVanderBoundary;

    vector<MatrixXd> GkVanderBoundary_comp(localSpace.Dimension);
    for (unsigned int d = 0; d < localSpace.Dimension; d++)
        GkVanderBoundary_comp[d] =
            GkVanderBoundary.block(0, d * numBoundaryQuadraturePoints, localSpace.Dimension * localSpace.Nk, numBoundaryQuadraturePoints);

    localSpace.GkVanderBoundaryTimesNormal =
        MatrixXd::Zero(localSpace.NumBoundaryBasisFunctions, localSpace.Dimension * localSpace.Nk);

    W2 = MatrixXd::Zero(localSpace.Nk, localSpace.NumBasisFunctions);

    B2Nabla = MatrixXd::Zero(localSpace.NkNabla, localSpace.NumBasisFunctions);

    const unsigned int nk2D = reference_element_data.NumDofs2D;

    unsigned int offsetQuadraturePoints = 0;

    const MatrixXd id = MatrixXd::Identity(nk2D, nk2D);

    localSpace.VanderBasisFunctionValuesOnFace.resize(numFaces);

    localSpace.FacesVanderInternal.resize(numFaces);
    for (unsigned int f = 0; f < numFaces; f++)
    {
        const double localDirection = facesNormalDirections[f] ? 1.0 : -1.0;
        const double globalDirection = faceNormalGlobalDirections[f] ? 1.0 : -1.0;

        const unsigned int numFaceQuadraturePoints = facesQuadrature[f].Weights.size();
        localSpace.FacesVanderInternal[f] =
            monomials2D.Vander(reference_element_data.Monomials2D, facesQuadrature[f].Points, facesCentroids[f], facesDiameters[f]);

        const VectorXd sqrtFacesQuadratureWeights = facesQuadrature[f].Weights.array().sqrt();

        const MatrixXd temp = sqrtFacesQuadratureWeights.asDiagonal() * localSpace.FacesVanderInternal[f];
        const MatrixXd faceHmatrix = temp.transpose() * temp;

        localSpace.VanderBasisFunctionValuesOnFace[f] =
            facesAreas[f] * localSpace.FacesVanderInternal[f] * faceHmatrix.llt().solve(id);

        W2.block(0, f * nk2D, localSpace.Nk, nk2D) =
            globalDirection *
            localSpace.VanderBoundary.block(offsetQuadraturePoints, 0, numFaceQuadraturePoints, localSpace.Nk).transpose() *
            facesQuadrature[f].Weights.asDiagonal() * localSpace.VanderBasisFunctionValuesOnFace[f];

        B2Nabla.block(0, f * nk2D, localSpace.NkNabla, nk2D) =
            globalDirection *
            localSpace.VanderBoundaryKp1.rightCols(localSpace.NkNabla)
                .transpose()
                .block(0, offsetQuadraturePoints, localSpace.NkNabla, numFaceQuadraturePoints) *
            facesQuadrature[f].Weights.asDiagonal() * localSpace.VanderBasisFunctionValuesOnFace[f];

        localSpace.GkVanderBoundaryTimesNormal.block(f * nk2D, 0, nk2D, 3 * localSpace.Nk) =
            (1.0 / facesAreas[f]) * localDirection * globalDirection * localSpace.FacesVanderInternal[f].transpose() *
            boundaryQuadratureWeights.segment(offsetQuadraturePoints, numFaceQuadraturePoints).asDiagonal() *
            (facesNormals[f](0) * GkVanderBoundary_comp[0].block(0, offsetQuadraturePoints, 3 * localSpace.Nk, numFaceQuadraturePoints) +
             facesNormals[f](1) * GkVanderBoundary_comp[1].block(0, offsetQuadraturePoints, 3 * localSpace.Nk, numFaceQuadraturePoints) +
             facesNormals[f](2) * GkVanderBoundary_comp[2].block(0, offsetQuadraturePoints, 3 * localSpace.Nk, numFaceQuadraturePoints))
                .transpose();

        offsetQuadraturePoints += numFaceQuadraturePoints;
    }
}
//****************************************************************************
} // namespace MCC
} // namespace VEM
} // namespace Polydim
