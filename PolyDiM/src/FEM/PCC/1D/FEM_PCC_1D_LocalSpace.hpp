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

#ifndef __FEM_PCC_1D_LocalSpace_HPP
#define __FEM_PCC_1D_LocalSpace_HPP

#include "FEM_PCC_1D_ReferenceElement.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_PCC_1D_Segment_Geometry final
{
    double Tolerance1D;

    Eigen::Vector3d Origin;
    Eigen::Vector3d Tangent;
    double Length;
};

struct FEM_PCC_1D_LocalSpace_Data final
{
    struct SegmentMapData
    {
        Eigen::Vector3d Origin;
        Eigen::Vector3d Tangent;
        double Length;
        double SquaredLength;
    };

    Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data::SegmentMapData MapData;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    Eigen::MatrixXd Dofs;
    std::vector<unsigned int> DofsMeshOrder;
    std::array<unsigned int, 3> Dof0DsIndex;
    std::array<unsigned int, 2> Dof1DsIndex;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
};

class FEM_PCC_1D_LocalSpace final
{
  private:
    Eigen::MatrixXd F(const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData, const Eigen::MatrixXd &x) const
    {
        Eigen::MatrixXd points(3, x.cols());

        for (unsigned int p = 0; p < x.cols(); ++p)
            points.col(p) << mapData.Origin + mapData.Tangent * x(0, p);

        return points;
    }

    Eigen::MatrixXd FInv(const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData, const Eigen::MatrixXd &x) const
    {
        Eigen::MatrixXd points(3, x.cols());

        for (unsigned int p = 0; p < x.cols(); ++p)
            points.col(p) << (x.col(p) - mapData.Origin).dot(mapData.Tangent) / mapData.SquaredLength;

        return points;
    }

    inline Eigen::VectorXd DetJ(const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData,
                                const Eigen::MatrixXd &x) const
    {
        return Eigen::VectorXd::Constant(x.cols(), mapData.Length);
    }

    Eigen::MatrixXd MapValues(const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space,
                              const Eigen::MatrixXd &referenceValues) const;

    std::vector<Eigen::MatrixXd> MapDerivativeValues(const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space,
                                                     const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData) const;

  public:
    Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data CreateLocalSpace(const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                                   const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry &segment) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                       const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                       const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space,
                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = FInv(local_space.MapData, points);

        Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Data &local_space,
        const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = FInv(local_space.MapData, points);

        Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
