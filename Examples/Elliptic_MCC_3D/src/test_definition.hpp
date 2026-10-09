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

#include <typeindex>
#include <unordered_map>

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_3D
{
namespace test
{
// ***************************************************************************
enum struct Test_Types
{
    Patch_Test = 1,
    Poisson_Polynomial_Problem = 2
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain() const = 0;
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

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain;

        domain.volume = 1.0;

        domain.vertices.resize(3, 8);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0;
        domain.vertices.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        // create edges
        domain.edges.setZero(2, 12);
        domain.edges.col(0) << 0, 1;
        domain.edges.col(1) << 1, 2;
        domain.edges.col(2) << 2, 3;
        domain.edges.col(3) << 3, 0;
        domain.edges.col(4) << 4, 5;
        domain.edges.col(5) << 5, 6;
        domain.edges.col(6) << 6, 7;
        domain.edges.col(7) << 7, 4;
        domain.edges.col(8) << 0, 4;
        domain.edges.col(9) << 1, 5;
        domain.edges.col(10) << 2, 6;
        domain.edges.col(11) << 3, 7;

        // create faces
        domain.faces.resize(6, Eigen::MatrixXi::Zero(2, 4));
        domain.faces[0].row(0) << 0, 1, 2, 3;
        domain.faces[0].row(1) << 0, 1, 2, 3;

        domain.faces[1].row(0) << 4, 5, 6, 7;
        domain.faces[1].row(1) << 4, 5, 6, 7;

        domain.faces[2].row(0) << 0, 3, 7, 4;
        domain.faces[2].row(1) << 3, 11, 7, 8;

        domain.faces[3].row(0) << 1, 2, 6, 5;
        domain.faces[3].row(1) << 1, 10, 5, 9;

        domain.faces[4].row(0) << 0, 1, 5, 4;
        domain.faces[4].row(1) << 0, 9, 4, 8;

        domain.faces[5].row(0) << 3, 2, 6, 7;
        domain.faces[5].row(1) << 2, 10, 6, 11;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types::Parallelepiped;

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {17, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {18, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {19, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {20, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {22, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {23, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {24, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {25, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}}};
    }

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -2.0)};
    }

    std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 16.0 / 3.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Constant(points.cols(), 4.0 / 3.0)};
    }

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(1).array() * points.row(2).array();
    }

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Constant(points.cols(), -2.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), -2.0),
                Eigen::VectorXd::Constant(points.cols(), 7.0)};
    };

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 17.0 / 3.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Constant(points.cols(), 5.0 / 3.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Constant(points.cols(), 2.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 5.0 / 3.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 2.0 / 3.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

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

        return -3.0 * second_derivatives + points.row(1).array().transpose() * points.row(0).array().transpose() *
                                               points.row(2).array().transpose() * solution;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

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
            return 4.0 * derivatives + 2.0 * solution;
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        if (order > 0)
        {
            derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for (int i = 0; i < order - 1; i++)
                derivatives = derivatives * polynomial;

            solution = derivatives * polynomial;
            derivatives *= order;
        }

        return {derivatives + solution, solution, -4.0 * derivatives - 2.0 * solution};
    }
};
// ***************************************************************************
struct Poisson_Polynomial_Problem final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain;

        domain.volume = 1.0;

        domain.vertices.resize(3, 8);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0;
        domain.vertices.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        // create edges
        domain.edges.setZero(2, 12);
        domain.edges.col(0) << 0, 1;
        domain.edges.col(1) << 1, 2;
        domain.edges.col(2) << 2, 3;
        domain.edges.col(3) << 3, 0;
        domain.edges.col(4) << 4, 5;
        domain.edges.col(5) << 5, 6;
        domain.edges.col(6) << 6, 7;
        domain.edges.col(7) << 7, 4;
        domain.edges.col(8) << 0, 4;
        domain.edges.col(9) << 1, 5;
        domain.edges.col(10) << 2, 6;
        domain.edges.col(11) << 3, 7;

        // create faces
        domain.faces.resize(6, Eigen::MatrixXi::Zero(2, 4));
        domain.faces[0].row(0) << 0, 1, 2, 3;
        domain.faces[0].row(1) << 0, 1, 2, 3;

        domain.faces[1].row(0) << 4, 5, 6, 7;
        domain.faces[1].row(1) << 4, 5, 6, 7;

        domain.faces[2].row(0) << 0, 3, 7, 4;
        domain.faces[2].row(1) << 3, 11, 7, 8;

        domain.faces[3].row(0) << 1, 2, 6, 5;
        domain.faces[3].row(1) << 1, 10, 5, 9;

        domain.faces[4].row(0) << 0, 1, 5, 4;
        domain.faces[4].row(1) << 0, 9, 4, 8;

        domain.faces[5].row(0) << 3, 2, 6, 7;
        domain.faces[5].row(1) << 2, 10, 6, 11;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types::Parallelepiped;

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {17, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {18, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {19, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {20, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {22, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {23, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {24, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {25, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}}};
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
                Eigen::VectorXd::Constant(points.cols(), 1.0)};
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
                Eigen::VectorXd::Constant(points.cols(), 1.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        return 128.0 *
               (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(2).array() * (1.0 - points.row(2).array()) +
                points.row(0).array() * (1.0 - points.row(0).array()) * points.row(2).array() * (1.0 - points.row(2).array()) +
                points.row(0).array() * (1.0 - points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()));
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

        return 64.0 * points.row(2).array() * (1.0 - points.row(2).array()) * points.row(1).array() *
                   (1.0 - points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) +
               1.7;
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return 64.0 * points.row(2).array() * (1.0 - points.row(2).array()) * points.row(1).array() *
                   (1.0 - points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) +
               1.7;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {-64.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()) *
                    points.row(2).array() * (1.0 - points.row(2).array()),
                -64.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) *
                    points.row(2).array() * (1.0 - points.row(2).array()),
                -64.0 * (1.0 - 2.0 * points.row(2).array()) * points.row(0).array() * (1.0 - points.row(0).array()) *
                    points.row(1).array() * (1.0 - points.row(1).array())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elliptic_MCC_3D
} // namespace examples
} // namespace Polydim

#endif
