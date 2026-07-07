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

/// @brief Abstract interface for a 2D primal conforming ZFEM (Zipped FEM) local space
class I_ZFEM_PCC_2D_LocalSpace
{
  public:
    /// \brief Class destructor.
    virtual ~I_ZFEM_PCC_2D_LocalSpace()
    {
    }

    /// @brief Build the ZFEM local space on a given polygon.
    ///
    /// @param reference_element_data Reference-element data (order, reference quantities).
    /// @param polygon                Geometry of the polygonal cell.
    /// @return The computed local-space data for the cell.
    virtual ZFEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const ZFEM_PCC_2D_Polygon_Geometry &polygon) const = 0;

    /// @brief Evaluate the basis functions at the cell's internal quadrature points.
    ///
    /// @param localSpace Local-space data for the cell (from CreateLocalSpace()).
    /// @return A (quadrature points \f$\times\f$ local DOFs) matrix of basis-function values.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    /// @brief Evaluate the basis-function derivatives at the cell's internal quadrature points.
    ///
    /// @param localSpace Local-space data for the cell (from CreateLocalSpace()).
    /// @return One (quadrature points \f$\times\f$ local DOFs) matrix per spatial direction.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ZFEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    /// @brief Evaluate the basis functions on a given edge.
    ///
    /// @param reference_element_data     Reference-element data.
    /// @param localSpace                 Local-space data for the cell (from CreateLocalSpace()).
    /// @param pointsCurvilinearCoordinates Evaluation points on the edge, in curvilinear coordinates.
    /// @return A matrix of edge basis-function values at the requested points.
    virtual Eigen::MatrixXd ComputeValuesOnEdge(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const ZFEM_PCC_2D_LocalSpace_Data &localSpace,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;
};
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
