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

#include <numbers>

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{
namespace test
{
// ***************************************************************************
enum struct Test_Types
{
    Patch_Test = 1,
    Poisson_Problem = 2 /// Test 1 in S. Berrone, S. Scialò, G. Teora. "The mixed virtual element discretization for
                        /// highly-anisotropic problems: the role of the boundary degrees of freedom". Mathematics in
                        /// Engineering, 2023, 5(6): 1-32. doi: 10.3934/mine.2023099
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static unsigned int order;

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
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}}};
    }

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Constant(points.cols(), -0.2),
                Eigen::VectorXd::Zero(points.cols())};
    }

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(1).array();
    }

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 2.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.6),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if (order > 1)
        {
            second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for (int i = 0; i < order - 2; i++)
                second_derivatives = second_derivatives * polynomial;

            solution = second_derivatives * polynomial * polynomial;
            second_derivatives *= order * (order - 1);
        }
        else if (order == 1)
            solution = polynomial;

        return -3.0 * second_derivatives + points.row(1).array().transpose() * points.row(0).array().transpose() * solution;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if (order > 0)
        {
            derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for (int i = 0; i < order - 1; i++)
                derivatives = derivatives * polynomial;

            solution = derivatives * polynomial;
            derivatives *= order;
        }

        switch (marker)
        {
        case 1:
            return 2.0 * derivatives + solution;
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if (order > 0)
        {
            derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for (int i = 0; i < order - 1; i++)
                derivatives = derivatives * polynomial;

            solution = derivatives * polynomial;
            derivatives *= order;
        }

        return {-derivatives + solution, -2.0 * derivatives - solution, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Poisson_Problem final : public I_Test
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
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}}};
    }

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Zero(points.cols());
    }

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        return 2.0 * std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) *
               sin(std::numbers::pi * points.row(1).array());
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        return sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array());
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {-std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                -std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elliptic_MCC_2D
} // namespace examples
} // namespace Polydim

#endif
