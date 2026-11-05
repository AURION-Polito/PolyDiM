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
namespace Parabolic_PCC_3D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Space_Test = 2,
    Time_Test = 3,
    Parabolic_Problem = 4
};

struct I_Test
{
    virtual ~I_Test()
    {
    }

    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_3D domain() const = 0;
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
    unsigned int space_order;
    unsigned int time_order;

    Patch_Test(const unsigned int space_order, const unsigned int time_order)
        : space_order(space_order), time_order(time_order)
    {
    }

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_3D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_3D domain;

        domain.spatial_domain.volume = 1.0;

        domain.spatial_domain.vertices.resize(3, 8);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0;
        domain.spatial_domain.vertices.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        // create edges
        domain.spatial_domain.edges.setZero(2, 12);
        domain.spatial_domain.edges.col(0) << 0, 1;
        domain.spatial_domain.edges.col(1) << 1, 2;
        domain.spatial_domain.edges.col(2) << 2, 3;
        domain.spatial_domain.edges.col(3) << 3, 0;
        domain.spatial_domain.edges.col(4) << 4, 5;
        domain.spatial_domain.edges.col(5) << 5, 6;
        domain.spatial_domain.edges.col(6) << 6, 7;
        domain.spatial_domain.edges.col(7) << 7, 4;
        domain.spatial_domain.edges.col(8) << 0, 4;
        domain.spatial_domain.edges.col(9) << 1, 5;
        domain.spatial_domain.edges.col(10) << 2, 6;
        domain.spatial_domain.edges.col(11) << 3, 7;

        // cre.spatial_domainate faces
        domain.spatial_domain.faces.resize(6, Eigen::MatrixXi::Zero(2, 4));
        domain.spatial_domain.faces[0].row(0) << 0, 1, 2, 3;
        domain.spatial_domain.faces[0].row(1) << 0, 1, 2, 3;

        domain.spatial_domain.faces[1].row(0) << 4, 5, 6, 7;
        domain.spatial_domain.faces[1].row(1) << 4, 5, 6, 7;

        domain.spatial_domain.faces[2].row(0) << 0, 3, 7, 4;
        domain.spatial_domain.faces[2].row(1) << 3, 11, 7, 8;

        domain.spatial_domain.faces[3].row(0) << 1, 2, 6, 5;
        domain.spatial_domain.faces[3].row(1) << 1, 10, 5, 9;

        domain.spatial_domain.faces[4].row(0) << 0, 1, 5, 4;
        domain.spatial_domain.faces[4].row(1) << 0, 9, 4, 8;

        domain.spatial_domain.faces[5].row(0) << 3, 2, 6, 7;
        domain.spatial_domain.faces[5].row(1) << 2, 10, 6, 11;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types::Parallelepiped;

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
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                //{11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
          {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {17, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {18, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {19, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {20, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                //{21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
          {21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {22, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {23, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {24, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {25, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                //{26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}}};
          {26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 2.0);
    }

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd source_space_term = Eigen::VectorXd::Constant(points.cols(), 6.0 * space_order * (space_order - 1));
        const Eigen::ArrayXd space_polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_space_order = space_order - 2;
        for (int i = 0; i < max_space_order; ++i)
            source_space_term.array() *= space_polynomial;

        Eigen::VectorXd source_time_term = Eigen::VectorXd::Constant(points.cols(), time_order);
        const Eigen::ArrayXd time_polynomial = Eigen::VectorXd::Constant(points.cols(), time_value);
        const int max_time_order = time_order - 1;
        for (int i = 0; i < max_time_order; ++i)
            source_time_term.array() *= time_polynomial;

        return source_time_term - source_space_term;
    }

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        return exact_solution(points, 0.0);
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    }

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {

        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), 2.0 * space_order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = space_order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        std::array<Eigen::VectorXd, 3> der = {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};

        switch (marker)
        {
        case 2: // co-normal derivatives on the left
            return - derivatives.array();
        case 4: // co-normal derivatives on the right
            return  derivatives.array();
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        const Eigen::ArrayXd space_polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;
        Eigen::VectorXd space_exact = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < space_order; ++i)
            space_exact.array() *= space_polynomial;

        const Eigen::ArrayXd time_polynomial = Eigen::VectorXd::Constant(points.cols(), time_value);
        Eigen::VectorXd time_exact = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < time_order; ++i)
            time_exact.array() *= time_polynomial;

        return space_exact + time_exact;
    }

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), space_order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = space_order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Parabolic_Problem final : public I_Test
{
    double K;

    Parabolic_Problem(const double &K) : K(K)
    {
    }

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_3D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_3D domain;

        domain.spatial_domain.volume = 1.0;

        domain.spatial_domain.vertices.resize(3, 8);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0;
        domain.spatial_domain.vertices.row(2) << 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        // create edges
        domain.spatial_domain.edges.setZero(2, 12);
        domain.spatial_domain.edges.col(0) << 0, 1;
        domain.spatial_domain.edges.col(1) << 1, 2;
        domain.spatial_domain.edges.col(2) << 2, 3;
        domain.spatial_domain.edges.col(3) << 3, 0;
        domain.spatial_domain.edges.col(4) << 4, 5;
        domain.spatial_domain.edges.col(5) << 5, 6;
        domain.spatial_domain.edges.col(6) << 6, 7;
        domain.spatial_domain.edges.col(7) << 7, 4;
        domain.spatial_domain.edges.col(8) << 0, 4;
        domain.spatial_domain.edges.col(9) << 1, 5;
        domain.spatial_domain.edges.col(10) << 2, 6;
        domain.spatial_domain.edges.col(11) << 3, 7;

        // cre.spatial_domainate faces
        domain.spatial_domain.faces.resize(6, Eigen::MatrixXi::Zero(2, 4));
        domain.spatial_domain.faces[0].row(0) << 0, 1, 2, 3;
        domain.spatial_domain.faces[0].row(1) << 0, 1, 2, 3;
        domain.spatial_domain.faces[1].row(0) << 4, 5, 6, 7;
        domain.spatial_domain.faces[1].row(1) << 4, 5, 6, 7;

        domain.spatial_domain.faces[2].row(0) << 0, 3, 7, 4;
        domain.spatial_domain.faces[2].row(1) << 3, 11, 7, 8;

        domain.spatial_domain.faces[3].row(0) << 1, 2, 6, 5;
        domain.spatial_domain.faces[3].row(1) << 1, 10, 5, 9;

        domain.spatial_domain.faces[4].row(0) << 0, 1, 5, 4;
        domain.spatial_domain.faces[4].row(1) << 0, 9, 4, 8;

        domain.spatial_domain.faces[5].row(0) << 3, 2, 6, 7;
        domain.spatial_domain.faces[5].row(1) << 2, 10, 6, 11;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types::Parallelepiped;

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
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {17, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {18, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {19, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {20, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {22, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {23, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {24, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {25, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), K);
    }

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points, const double &time_value) const
    {
        return (1.0 + K * 2.0 * std::numbers::pi * std::numbers::pi) * sin(std::numbers::pi * points.row(0).array()) *
               sin(std::numbers::pi * points.row(1).array()) * exp(time_value);
    }

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        return exact_solution(points, 0.0);
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    }

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
        return sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()) * exp(time_value) +
               (K - 1.0);
    }

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
} // namespace Parabolic_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
