#ifndef __VEM_MCC_3D_Pressure_LocalSpace_HPP
#define __VEM_MCC_3D_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_3D_Pressure_LocalSpace.hpp"
#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_MCC_3D_LocalSpace_Data.hpp"
#include "VEM_MCC_Utilities.hpp"
#include "VEM_Monomials_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{

class VEM_MCC_3D_Pressure_LocalSpace final : public I_VEM_MCC_3D_Pressure_LocalSpace
{
  private:
    VEM_MCC_Utilities<3> utilities;
    Utilities::VEM_Monomials_3D monomials;

    void InitializeProjectorsComputation(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const;

  public:
    VEM_MCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                         const VEM_MCC_3D_Polyhedron_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                       const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                       const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                    const VEM_MCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
