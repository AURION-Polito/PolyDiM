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

#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Elliptic_Polynomial_Problem = 2, /// Test 1: S. Berrone, G. Teora, F. Vicini, "Improving high-order VEM stability on
    /// badly-shaped elements", doi: https://doi.org/10.1016/j.matcom.2023.10.003.
    Elliptic_Problem = 3
};

struct PDE_Time_Domain_2D
{
    std::array<double, 2> time_domain;
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D spatial_domain;
};

struct I_Test
{
    virtual PDE_Time_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                      const Eigen::MatrixXd &points,
                                                      const double &time_value) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                    const Eigen::MatrixXd &points,
                                                    const double &time_value) const = 0;
    virtual Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points,
                                                                     const double &time_value) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static unsigned int order;

    PDE_Time_Domain_2D domain() const
    {
      PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = { 0.0, 1.0 };

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        source_term *= -1.0;
        source_term.array() += 1.0;

        return source_term;
    };

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        return exact_solution(points, 0.0);
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {

        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        std::array<Eigen::VectorXd, 3> der = {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};

        switch (marker)
        {
        case 2:
            return order * derivatives.array();
        case 4:
            return order * derivatives.array();
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        result.array() += time_value;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Elliptic_Polynomial_Problem final : public I_Test
{
    PDE_Time_Domain_2D domain() const
    {
        PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = { 0.0, 1.0 };

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

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        const double k = 1.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        return exact_solution(points, 0.0);
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        switch (marker)
        {
        case 2: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 4: // co-normal derivatives on the left
            return -16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return {16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Elliptic_Problem final : public I_Test
{
    PDE_Time_Domain_2D domain() const
    {
        PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = { 0.0, 1.0 };

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

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        const double k = 2.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return 16.0 * std::numbers::pi * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
               sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        return exact_solution(points, 0.0);
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return sin(2.0 * std::numbers::pi * points.row(0).array()) * sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return {2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * points.row(0).array()) *
                    sin(2.0 * std::numbers::pi * points.row(1).array()),
                2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
                    cos(2.0 * std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
