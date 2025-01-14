#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

#include <typeindex>
#include <unordered_map>

namespace Polydim
{
namespace examples
{
namespace Elastic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Elasticity_Benchmark_1 = 2 /// Test 5.2a Da Veiga - Lovadina 2015
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 2> lame_coefficients(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker,
                                                                     const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker,
                                                                   const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_displacement(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 9> exact_derivatives_displacement(const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static int order;

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        if (order < 2)
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 2; i++)
            result = result * polynomial;

        return {-6.0 * order * (order - 1) * result, -6.0 * order * (order - 1) * result, Eigen::VectorXd::Zero(points.cols())};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return {result, result, Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
            result = result * polynomial;

        switch (marker)
        {
        case 2: // co-normal derivatives on the bottom
            return {order * result, order * result, Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    std::array<Eigen::VectorXd, 2> lame_coefficients(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        return {Eigen::ArrayXd::Constant(points.cols(), 1.0), Eigen::ArrayXd::Constant(points.cols(), 1.0)};
    }

    std::array<Eigen::VectorXd, 3> exact_displacement(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return {result, result, Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_displacement(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
            result = result * polynomial;

        result *= order;

        return {result,
                result,
                Eigen::VectorXd::Zero(points.cols()),
                result,
                result,
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
