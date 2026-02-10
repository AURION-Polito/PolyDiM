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

#ifndef __I_ZFEM_PCC_2D_LocalSpace_HPP
#define __I_ZFEM_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_ZFEM_PCC_2D_ReferenceElement.hpp"
#include "ZFEM_PCC_2D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{

class I_ZFEM_PCC_2D_LocalSpace
{
  public:
    /// \brief Class destructor.
    virtual ~I_ZFEM_PCC_2D_LocalSpace()
    {
    }

    virtual ZFEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const ZFEM_PCC_2D_Polygon_Geometry &polygon) const = 0;

    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    virtual Eigen::MatrixXd ComputeValuesOnEdge(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const ZFEM_PCC_2D_LocalSpace_Data &localSpace,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;
};
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
