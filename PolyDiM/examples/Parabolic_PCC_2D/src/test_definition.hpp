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
    Space_Test = 2,
    Time_Test = 3,
    Parabolic_Problem = 4 /// Test 1: G. Vacca, L. Beirão da Veiga, "Virtual element methods for parabolic problems on
                          /// polygonal meshes", doi: https://doi.org/10.1002/num.21982.
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const = 0;
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
    static unsigned int space_order;
    static unsigned int time_order;

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = {0.0, 1.0};

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
        Eigen::VectorXd source_space_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * space_order * (space_order - 1));
        const Eigen::ArrayXd space_polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_space_order = space_order - 2;
        for (int i = 0; i < max_space_order; ++i)
            source_space_term.array() *= space_polynomial;

        Eigen::VectorXd source_time_term = Eigen::VectorXd::Constant(points.cols(), time_order);
        const Eigen::ArrayXd time_polynomial = Eigen::VectorXd::Constant(points.cols(), time_value);
        const int max_time_order = time_order - 1;
        for (int i = 0; i < max_time_order; ++i)
            source_time_term.array() *= time_polynomial;

        return source_time_term - source_space_term;
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

        const int max_order = space_order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        std::array<Eigen::VectorXd, 3> der = {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};

        switch (marker)
        {
        case 2:
            return space_order * derivatives.array();
        case 4:
            return space_order * derivatives.array();
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        const Eigen::ArrayXd space_polynomial = points.row(0).array() + points.row(1).array() + 0.5;
        Eigen::VectorXd space_exact = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < space_order; ++i)
            space_exact.array() *= space_polynomial;

        const Eigen::ArrayXd time_polynomial = Eigen::VectorXd::Constant(points.cols(), time_value);
        Eigen::VectorXd time_exact = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < time_order; ++i)
            time_exact.array() *= time_polynomial;

        return space_exact + time_exact;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), space_order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = space_order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Parabolic_Problem final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = {0.0, 1.0};

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
        return (1.0 + 2.0 * std::numbers::pi * std::numbers::pi) * sin(std::numbers::pi * points.row(0).array()) *
               sin(std::numbers::pi * points.row(1).array()) * exp(time_value);
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
        return sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()) * exp(time_value);
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return {std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) *
                    sin(std::numbers::pi * points.row(1).array()) * exp(time_value),
                std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) *
                    cos(std::numbers::pi * points.row(1).array()) * exp(time_value),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
