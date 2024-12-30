#ifndef __FEM_PCC_1D_LocalSpace_HPP
#define __FEM_PCC_1D_LocalSpace_HPP

#include "FEM_PCC_1D_ReferenceElement.hpp"
#include "MapTriangle.hpp"

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

    SegmentMapData MapData;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    Eigen::MatrixXd Dofs;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
};

/// \brief Interface used to FEM Values computation
class FEM_PCC_1D_LocalSpace final
{
  private:
    Eigen::MatrixXd F(const FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData, const Eigen::MatrixXd &x) const
    {
        Eigen::MatrixXd points(3, x.cols());

        for (unsigned int p = 0; p < x.cols(); ++p)
            points.col(p) << mapData.Origin + mapData.Tangent * x(0, p);

        return points;
    }

    Eigen::MatrixXd FInv(const FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData, const Eigen::MatrixXd &x) const
    {
        Eigen::MatrixXd points(3, x.cols());

        for (unsigned int p = 0; p < x.cols(); ++p)
            points.col(p) << (x.col(p) - mapData.Origin).dot(mapData.Tangent) / mapData.SquaredLength;

        return points;
    }

    inline Eigen::VectorXd DetJ(const FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData, const Eigen::MatrixXd &x) const
    {
        return Eigen::VectorXd::Constant(x.cols(), mapData.Length);
    }

    /// \brief map basis function values on element with correct order
    inline Eigen::MatrixXd MapValues(const FEM_PCC_1D_LocalSpace_Data &, const Eigen::MatrixXd &referenceValues) const
    {
        return referenceValues;
    }

    /// \brief map basis function derivative values on element with correct order
    inline std::vector<Eigen::MatrixXd> MapDerivativeValues(const FEM_PCC_1D_LocalSpace_Data &local_space,
                                                            const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
    {
        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(1);

        basisFunctionsDerivativeValues[0] = referenceDerivateValues[0] / local_space.MapData.Length;

        return basisFunctionsDerivativeValues;
    }

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const FEM_PCC_1D_LocalSpace_Data::SegmentMapData &mapData) const;

  public:
    FEM_PCC_1D_LocalSpace_Data CreateLocalSpace(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                const FEM_PCC_1D_Segment_Geometry &segment) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_PCC_1D_LocalSpace_Data &local_space) const
    {
        return MapValues(local_space, reference_element_data.ReferenceBasisFunctionValues);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_PCC_1D_LocalSpace_Data &local_space) const
    {
        return MapDerivativeValues(local_space, reference_element_data.ReferenceBasisFunctionDerivativeValues);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                       const FEM_PCC_1D_LocalSpace_Data &local_space,
                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = FInv(local_space.MapData, points);

        FEM_PCC_1D_ReferenceElement reference_element;

        return MapValues(local_space, reference_element.EvaluateBasisFunctions(referencePoints, reference_element_data));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                                              const FEM_PCC_1D_LocalSpace_Data &local_space,
                                                                              const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints = FInv(local_space.MapData, points);

        FEM_PCC_1D_ReferenceElement reference_element;

        return MapDerivativeValues(local_space, reference_element.EvaluateBasisFunctionDerivatives(referencePoints, reference_element_data));
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
