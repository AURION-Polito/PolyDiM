// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __VEM_PCC_3D_LocalSpace_HPP
#define __VEM_PCC_3D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_LocalSpace.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{
class VEM_PCC_3D_LocalSpace final : public I_VEM_PCC_3D_LocalSpace
{
private:
    VEM_PCC_Utilities<3> utilities;
    Utilities::VEM_Monomials_3D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polyhedronVertices,
                                         const Eigen::MatrixXi &polyhedronEdges,
                                         const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronMeasure,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &edgeInternalQuadraturePoints,
                                         VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                        const double &polyhedronMeasure,
                        const double &polyhedronDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const Eigen::VectorXd &boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                        VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                          const double &polyhedronMeasure,
                                          const double &polyhedronDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polyhedronMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        utilities.ComputeL2Projectors(polyhedronMeasure,
                                      localSpace.Order,
                                      localSpace.Nkm1,
                                      localSpace.NumProjectorBasisFunctions,
                                      localSpace.NumInternalBasisFunctions,
                                      localSpace.NumBasisFunctions,
                                      localSpace.Hmatrix,
                                      localSpace.PiNabla,
                                      localSpace.Cmatrix,
                                      localSpace.Pi0km1,
                                      localSpace.Pi0k);
    };

    void ComputeFaceProjectors(const VEM_PCC_2D_LocalSpace &faceVemValues,
                               const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                               const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                               const Eigen::VectorXd &boundaryQuadratureWeights,
                               VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputePolynomialsDofs(const double &polytopeMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const;

public:
    VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
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



    inline Eigen::MatrixXd ComputeDRecipeStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                      const ProjectionTypes &projectionType,
                                                      const Eigen::MatrixXd &coercivity_matrix,
                                                      const Eigen::VectorXd &vector_coefficients) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla:
            return utilities.ComputeDRecipeStabilizationMatrix(localSpace.PiNabla, coercivity_matrix, vector_coefficients, localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDRecipeStabilizationMatrix(localSpace.Pi0k, coercivity_matrix, vector_coefficients, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(projectionType,
                                                              localSpace.Nkm1,
                                                              localSpace.VanderInternalDerivatives,
                                                              localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, localSpace, points),
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                const ProjectionTypes &projectionType,
                                                                const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(
            projectionType,
            localSpace.Nkm1,
            localSpace.Pi0km1Der,
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  localSpace.Centroid,
                                                  localSpace.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            localSpace.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                             const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           localSpace.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                               const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        return utilities.ComputeValuesOnEdge(localSpace.EdgeInternalPoints.transpose(), reference_element_data.Order, localSpace.EdgeBasisCoefficients, pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
