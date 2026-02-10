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

#ifndef __ZFEM_PCC_2D_LocalSpace_HPP
#define __ZFEM_PCC_2D_LocalSpace_HPP

#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "I_ZFEM_PCC_2D_LocalSpace.hpp"
#include "Monomials_2D.hpp"
#include "ZFEM_PCC_Utilities.hpp"

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{

class ZFEM_PCC_2D_LocalSpace final : public I_ZFEM_PCC_2D_LocalSpace
{
  private:
    ZFEM_PCC_Utilities ZFEM_utilities;
    Utilities::Monomials_2D monomials;
    FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace fem_local_space;

    void PolygonFineNodes(const unsigned int num_vertices,
                          const unsigned int &NumBasisFunctions,
                          const unsigned int &NumTotalBasisFunctions,
                          const unsigned int &NumVirtualBasisFunctions,
                          const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
                          const Eigen::MatrixXi &local_to_total,
                          Eigen::MatrixXd &CoarseNodes,
                          Eigen::MatrixXd &VirtualNodes,
                          Eigen::MatrixXd &FinerNodes) const;

  public:
    ZFEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                 const ZFEM_PCC_2D_Polygon_Geometry &polygon) const;

    void ComputePolynomialsDofs(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                const Eigen::Vector3d &internal_points,
                                const double &diameter,
                                ZFEM_PCC_2D_LocalSpace_Data &localSpace) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return localSpace.ZFEM_basis_functions_values;
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return localSpace.ZFEM_basis_functions_derivative_values;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const ZFEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const std::vector<Eigen::MatrixXd> &points) const
    {
        const Eigen::MatrixXd total_fem_basis_functions_values =
            ZFEM_utilities.ComputeFEMBasisFunctionsValues(reference_element_data,
                                                          localSpace.NumBasisFunctions,
                                                          localSpace.NumVirtualBasisFunctions,
                                                          localSpace.fem_local_space_data,
                                                          localSpace.local_to_total,
                                                          points);

        return ZFEM_utilities.ComputeBasisFunctionsValues(localSpace.NumBasisFunctions,
                                                          localSpace.NumVirtualBasisFunctions,
                                                          total_fem_basis_functions_values,
                                                          localSpace.VirtualWeights);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                              const ZFEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const std::vector<Eigen::MatrixXd> &points) const
    {

        const std::vector<Eigen::MatrixXd> total_fem_basis_functions_derivative_values =
            ZFEM_utilities.ComputeFEMBasisFunctionsDerivativeValues(localSpace.Dimension,
                                                                    reference_element_data,
                                                                    localSpace.NumBasisFunctions,
                                                                    localSpace.NumVirtualBasisFunctions,
                                                                    localSpace.fem_local_space_data,
                                                                    localSpace.local_to_total,
                                                                    points);

        return ZFEM_utilities.ComputeBasisFunctionsDerivativeValues(localSpace.Dimension,
                                                                    localSpace.NumBasisFunctions,
                                                                    localSpace.NumVirtualBasisFunctions,
                                                                    total_fem_basis_functions_derivative_values,
                                                                    localSpace.VirtualWeights);
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const ZFEM_PCC_2D_LocalSpace_Data &localSpace,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {

        return ZFEM_utilities.ComputeValuesOnEdge(localSpace.ReferenceEdgeDOFsInternalPoints,
                                                  reference_element_data.Order,
                                                  localSpace.EdgeBasisCoefficients,
                                                  pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
