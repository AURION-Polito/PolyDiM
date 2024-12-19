#ifndef __VEM_MCC_3D_Velocity_LocalSpace_HPP
#define __VEM_MCC_3D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_MCC_Utilities.hpp"
#include "VEM_MCC_Velocity_LocalSpace_Data.hpp"
#include "VEM_Monomials_2D.hpp"
#include "VEM_Monomials_3D.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
/// \brief Class used for computing values of basis functions of 3D
/// Mixed Conforming Constant degree Virtual Element Methods.
class VEM_MCC_3D_Velocity_LocalSpace final : public I_VEM_MCC_3D_Velocity_LocalSpace
{
private:
    MCC::VEM_MCC_Utilities<3> utilities;
    Monomials::VEM_Monomials_3D monomials3D;
    Monomials::VEM_Monomials_2D monomials2D;

    void InitializeProjectorsComputation(const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const unsigned int &numFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         VEM_MCC_Velocity_LocalSpace_Data &localSpace) const;

    inline void ComputeStabilizationMatrix(const double &polyhedronMeasure,
                                           VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
    {
        localSpace.StabMatrix =
            utilities.ComputeStabilizationMatrix(localSpace.Pi0k, polyhedronMeasure, localSpace.Dmatrix);
    }

    void ComputeL2Projectors(const double &polyhedronMeasure,
                             const Eigen::VectorXd &internalQuadratureWeights,
                             const Eigen::MatrixXd &B2Nabla,
                             VEM_MCC_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       VEM_MCC_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                 const unsigned int &numFaces,
                                 const std::vector<Eigen::Vector3d> &facesNormals,
                                 const std::vector<bool> &facesNormalDirections,
                                 const std::vector<bool> &faceNormalGlobalDirections,
                                 const std::vector<Eigen::Vector3d> &facesCentroids,
                                 const std::vector<double> &facesAreas,
                                 const std::vector<double> &facesDiameters,
                                 const std::vector<Gedim::Quadrature::QuadratureData> &facesQuadrature,
                                 const Eigen::VectorXd &boundaryQuadratureWeights,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 VEM_MCC_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const double &polytopeMeasure, VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
    {
        localSpace.Dmatrix = utilities.ComputePolynomialBasisDofs(polytopeMeasure,
                                                                  localSpace.Order,
                                                                  localSpace.Nk,
                                                                  localSpace.NumBoundaryBasisFunctions,
                                                                  localSpace.NumNablaInternalBasisFunctions,
                                                                  localSpace.NumBigOPlusInternalBasisFunctions,
                                                                  localSpace.NumBasisFunctions,
                                                                  localSpace.GkVanderBoundaryTimesNormal,
                                                                  localSpace.Gmatrix);
    };

public:
    VEM_MCC_Velocity_LocalSpace_Data CreateLocalSpace(
        const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
        const VEM_MCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(
        const VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
    {
        const unsigned int numQuadrature = localSpace.InternalQuadrature.Points.cols();
        const Eigen::MatrixXd temp = localSpace.GkVanderInternal.transpose() * localSpace.Pi0k;
        std::vector<Eigen::MatrixXd> result(localSpace.Dimension,
                                            Eigen::MatrixXd::Zero(localSpace.Dimension, localSpace.NumBasisFunctions));

        for (unsigned int d = 0; d < localSpace.Dimension; d++)
            result[d] = temp.middleRows(numQuadrature * d, numQuadrature);

        return result;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(
        const VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Vmatrix;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(
        const VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
        const VEM_MCC_Velocity_LocalSpace_Data &localSpace,
        const VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
        const Eigen::MatrixXd &points) const
    {
        return monomials3D.Vander(reference_element_data.MonomialsKp1, points, polyhedron.Centroid, polyhedron.Diameter)
            .leftCols(localSpace.Nk);
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
