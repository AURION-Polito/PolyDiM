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

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_BulkFace_2D
{
namespace test
{
enum struct Test_Types
{
    Elliptic_Problem = 1, /// Experiment 1: M. Frittelli, A. Madzvamuse, I. Sgura, "Bulk-surface virtual element method
    /// for systems of PDEs in two-space dimensions", doi: https://doi.org/10.1007/s00211-020-01167-3.
    Parabolic_Problem = 2, /// Experiment 2: M. Frittelli, A. Madzvamuse, I. Sgura, "Bulk-surface virtual element method
    /// for systems of PDEs in two-space dimensions", doi: https://doi.org/10.1007/s00211-020-01167-3.
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const = 0;

    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_2D() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_1D() const = 0;

    virtual Eigen::VectorXd diffusion_term(const unsigned int dimension,
                                           const unsigned int id_domain,
                                           const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::VectorXd reaction_term(const unsigned int dimension,
                                          const unsigned int id_domain,
                                          const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::VectorXd source_term(const unsigned int dimension,
                                        const unsigned int id_domain,
                                        const Eigen::MatrixXd &points,
                                        const double &time_value) const = 0;

    virtual Eigen::VectorXd initial_solution(const unsigned int dimension,
                                             const unsigned int id_domain,
                                             const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::VectorXd exact_solution(const unsigned int dimension,
                                           const unsigned int id_domain,
                                           const Eigen::MatrixXd &points,
                                           const double &time_value) const = 0;

    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int dimension,
                                                    const unsigned int id_domain,
                                                    const unsigned int marker,
                                                    const Eigen::MatrixXd &points,
                                                    const double &time_value) const = 0;

    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int dimension,
                                                                     const unsigned int id_domain,
                                                                     const Eigen::MatrixXd &points,
                                                                     const double &time_value) const = 0;
    virtual Eigen::VectorXd alpha(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd beta(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Elliptic_Problem final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain;
        domain.time_domain = {};

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_2D() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_1D() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        double kmax = 0.0;
        switch (dimension)
        {
        case 2:
            kmax = 1.0;
            break;
        case 1:
            kmax = 1.0;
            break;
        default:
            throw std::runtime_error("not valid dimension");
        }

        return Eigen::VectorXd::Constant(points.cols(), kmax);
    };

    Eigen::VectorXd reaction_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        double kmax = 0.0;
        switch (dimension)
        {
        case 2:
            kmax = 1.0;
            break;
        case 1:
            kmax = 1.0;
            break;
        default:
            throw std::runtime_error("not valid dimension");
        }

        return Eigen::VectorXd::Constant(points.cols(), kmax);
    };

    Eigen::VectorXd source_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points, const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Constant(points.cols(), 0.0).array() +
                   exact_solution(dimension, id_domain, points, time_value).array();
        }
        case 1: {
            return 6.0 * points.row(0).array() * points.row(1).array() + 3.5 * points.row(0).array() * points.row(1).array();
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    };

    Eigen::VectorXd initial_solution(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return exact_solution(dimension, id_domain, points, 0.0);
    };

    Eigen::VectorXd exact_solution(const unsigned int dimension,
                                   const unsigned int id_domain,
                                   const Eigen::MatrixXd &points,
                                   const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return points.row(0).array() * points.row(1).array();
        }
        case 1: {
            return 1.5 * points.row(0).array() * points.row(1).array();
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int dimension,
                                            const unsigned int id_domain,
                                            const unsigned int marker,
                                            const Eigen::MatrixXd &points,
                                            const double &time_value) const
    {
        throw std::runtime_error("not valid condition");
    }

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int dimension,
                                                             const unsigned int id_domain,
                                                             const Eigen::MatrixXd &points,
                                                             const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return {points.row(1).array(), points.row(0).array(), Eigen::VectorXd::Constant(points.cols(), 0.0)};
        }
        case 1: {
            return {1.5 * points.row(1).array(), 1.5 * points.row(0).array(), Eigen::VectorXd::Constant(points.cols(), 0.0)};
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }

    Eigen::VectorXd alpha(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Constant(points.cols(), 1.0);
        }
        case 1: {
            return Eigen::VectorXd::Constant(points.cols(), 1.0);
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }

    Eigen::VectorXd beta(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Constant(points.cols(), 2.0);
        }
        case 1: {
            return Eigen::VectorXd::Constant(points.cols(), 2.0);
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }
};
// ***************************************************************************
struct Parabolic_Problem final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain;
        domain.time_domain = {0.0, 1.0};

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_2D() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info_1D() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        double kmax = 0.0;
        switch (dimension)
        {
        case 2:
            kmax = 1.0;
            break;
        case 1:
            kmax = 0.25;
            break;
        default:
            throw std::runtime_error("not valid dimension");
        }

        return Eigen::VectorXd::Constant(points.cols(), kmax);
    };

    Eigen::VectorXd reaction_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        double kmax = 0.0;
        switch (dimension)
        {
        case 2:
            kmax = 1.0;
            break;
        case 1:
            kmax = 0.0;
            break;
        default:
            throw std::runtime_error("not valid dimension");
        }

        return Eigen::VectorXd::Constant(points.cols(), kmax);
    };

    Eigen::VectorXd source_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points, const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Zero(points.cols());
        }
        case 1: {
            return Eigen::VectorXd::Zero(points.cols());
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    };

    Eigen::VectorXd initial_solution(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return exact_solution(dimension, id_domain, points, 0.0);
    };

    Eigen::VectorXd exact_solution(const unsigned int dimension,
                                   const unsigned int id_domain,
                                   const Eigen::MatrixXd &points,
                                   const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return points.row(0).array() * points.row(1).array() * exp(-time_value);
        }
        case 1: {
            return 1.5 * points.row(0).array() * points.row(1).array() * exp(-time_value);
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int dimension,
                                            const unsigned int id_domain,
                                            const unsigned int marker,
                                            const Eigen::MatrixXd &points,
                                            const double &time_value) const
    {
        throw std::runtime_error("not valid condition");
    }

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int dimension,
                                                             const unsigned int id_domain,
                                                             const Eigen::MatrixXd &points,
                                                             const double &time_value) const
    {
        switch (dimension)
        {
        case 2: {
            return {points.row(1).array() * exp(-time_value),
                    points.row(0).array() * exp(-time_value),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        }
        case 1: {
            return {1.5 * points.row(1).array() * exp(-time_value),
                    1.5 * points.row(0).array() * exp(-time_value),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }

    Eigen::VectorXd alpha(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Constant(points.cols(), 0.0);
        }
        case 1: {
            return Eigen::VectorXd::Constant(points.cols(), 2.0);
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }

    Eigen::VectorXd beta(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        switch (dimension)
        {
        case 2: {
            return Eigen::VectorXd::Constant(points.cols(), 4.0 / 3.0);
        }
        case 1: {
            return Eigen::VectorXd::Constant(points.cols(), 4.0 / 3.0);
        }
        default:
            throw std::runtime_error("not valid dimension");
        }
    }
};
// ***************************************************************************
} // namespace test
} // namespace Parabolic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim

#endif
