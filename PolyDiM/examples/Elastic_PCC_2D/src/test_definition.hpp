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
    LinearElasticity = 2
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
        const double mu = 1.0;
        const double lambda = 1.0;
        return {Eigen::ArrayXd::Constant(points.cols(), mu), Eigen::ArrayXd::Constant(points.cols(), lambda)};
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
struct LinearElasticity final : public I_Test
{
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
        const double mu = 1.5;
        const double lambda = 3.0;

        const Eigen::VectorXd u_1_xx = -32.0 * points.row(1).array() * (1.0 - points.row(1).array());
        const Eigen::VectorXd u_1_xy = 16.0 * (1.0 - 2.0 * points.row(1).array()) * (1.0 - 2.0 * points.row(0).array());
        const Eigen::VectorXd u_1_yy = -32.0 * points.row(0).array() * (1.0 - points.row(0).array());

        const Eigen::VectorXd u_2_xx = 5.0 * u_1_xx;
        const Eigen::VectorXd u_2_xy = 5.0 * u_1_xy;
        const Eigen::VectorXd u_2_yy = 5.0 * u_1_yy;

        return {-(2.0 * mu * (u_1_xx + 0.5 * (u_1_yy + u_2_xy)) + lambda * (u_1_xx + u_2_xy)),
                -(2.0 * mu * (0.5 * (u_1_xy + u_2_xx) + u_2_yy) + lambda * (u_1_xy + u_2_yy)),
                Eigen::VectorXd::Zero(points.cols())};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::VectorXd g = 16.0 * points.row(0).array() * (1.0 - points.row(0).array()) * points.row(1).array() *
                                      (1.0 - points.row(1).array()) +
                                  1.1;
        return {g, 5.0 * g, Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    std::array<Eigen::VectorXd, 2> lame_coefficients(const Eigen::MatrixXd &points) const
    {
        const double mu = 1.5;
        const double lambda = 3.0;
        return {Eigen::ArrayXd::Constant(points.cols(), mu), Eigen::ArrayXd::Constant(points.cols(), lambda)};
    }

    std::array<Eigen::VectorXd, 3> exact_displacement(const Eigen::MatrixXd &points) const
    {
        const Eigen::VectorXd g = 16.0 * points.row(0).array() * (1.0 - points.row(0).array()) * points.row(1).array() *
                                      (1.0 - points.row(1).array()) +
                                  1.1;

        return {g, 5.0 * g, Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_displacement(const Eigen::MatrixXd &points) const
    {
        return {16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols()),
                5.0 * 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                5.0 * 16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
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
