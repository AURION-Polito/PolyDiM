#ifndef __VEM_DF_PCC_3D_Velocity_LocalSpace_HPP
#define __VEM_DF_PCC_3D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_DF_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_3D_Velocity_LocalSpace final : public I_VEM_DF_PCC_3D_Velocity_LocalSpace
{
  private:
    VEM_DF_PCC_Utilities<3> utilities;
    Monomials::VEM_Monomials_3D monomials;
    Monomials::VEM_GBasis_3D g_basis;

    void InitializeProjectorsComputation(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polyhedronVertices,
                                         const Eigen::MatrixXi &polyhedronEdges,
                                         const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &boundaryQuadratureKLPoints,
                                         const Eigen::MatrixXd &edgeInternalQuadraturePoints,
                                         VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                       const double &polyhedronMeasure,
                                       const double &polyhedronDiameter,
                                       const std::vector<bool> &faceNormalGlobalDirections,
                                       const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                       VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const double &polyhedronMeasure,
                                    const double &polyhedronDiameter,
                                    const Eigen::VectorXd &internalQuadratureWeights,
                                    const Eigen::VectorXd &boundaryQuadratureWeights,
                                    VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeCMatrixkm2(const double &polyhedronMeasure,
                           const double &polyhedronDiameter,
                           const Eigen::VectorXd &boundaryQuadratureKLWeights,
                           VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                        const double &polyhedronMeasure,
                        const double &polyhedronDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                        VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights, VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                          const double &polyhedronDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeFaceProjectors(const PCC::VEM_PCC_2D_LocalSpace &faceVemValues,
                               const std::vector<PCC::VEM_PCC_2D_LocalSpace_Data> &facesLocalSpace,
                               const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                               const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                               const std::vector<std::array<Eigen::Vector3d, 2>> &facesTangents,
                               const std::vector<std::array<bool, 2>> &facesTangentsGlobalDirection,
                               const std::vector<Eigen::Vector3d> &facesNormals,
                               const std::vector<bool> &facesNormalDirections,
                               const std::vector<bool> &facesNormalGlobalDirections,
                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                               const Eigen::MatrixXd &boundaryQuadratureKLPoints,
                               VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const;

  public:
    VEM_DF_PCC_3D_Velocity_LocalSpace_Data CreateLocalSpace(const PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                            const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data_3D,
                                                            const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                            const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.PiNabla, localSpace.Diameter, localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.Measure, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                    const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType,
                                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, polyhedron, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                              const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, polyhedron, points),
            ComputePolynomialsDerivativeValues(reference_element_data, polyhedron, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  polyhedron.Centroid,
                                                  polyhedron.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polyhedron.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, polyhedron, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                             const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polyhedron.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, polyhedron, points));
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsDivergenceValues(localSpace.Nkm1, localSpace.VanderInternal, localSpace.Vmatrix);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
