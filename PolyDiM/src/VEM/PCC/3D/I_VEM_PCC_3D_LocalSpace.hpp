#ifndef __I_VEM_PCC_3D_LocalSpace_HPP
#define __I_VEM_PCC_3D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{

class I_VEM_PCC_3D_LocalSpace
{
  public:
    virtual ~I_VEM_PCC_3D_LocalSpace()
    {
    }

    virtual VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                        const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                        const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                        const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const = 0;

    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                 const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                        const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType,
                                                        const Eigen::MatrixXd &points) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                               const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType,
                                                                               const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                 const ProjectionTypes &projectionType,
                                                                 const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                     const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                            const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                              const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                              const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                const Eigen::VectorXd &edgeInternalPoints,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
