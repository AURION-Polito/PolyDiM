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
    Patch_Test = 1,
    Elliptic_Problem = 2, /// Experiment 1: M. Frittelli, A. Madzvamuse, I. Sgura, "Bulk-surface virtual element method
    /// for systems of PDEs in
    /// two-space dimensions", doi: https://doi.org/10.1007/s00211-020-01167-3.
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const = 0;

    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;

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

    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int dimension,
                                                                     const unsigned int id_domain,
                                                                     const Eigen::MatrixXd &points,
                                                                     const double &time_value) const = 0;
    virtual Eigen::VectorXd alpha(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd beta(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
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
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 0.0);
    };

    Eigen::VectorXd reaction_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    Eigen::VectorXd source_term(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 0.0);

        switch (dimension)
        {
        case 2:
            source_term = /*-source_term +*/ exact_solution(dimension, id_domain, points, time_value).array() + 1.0;
            break;
        case 1:
            source_term = /*-source_term +*/ exact_solution(dimension, id_domain, points, time_value).array() + 1.0;
            break;
        default:
            throw std::runtime_error("not valid dimension");
        }
        return source_term;
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
        const Eigen::ArrayXd polynomial = (points.row(0).array() + points.row(1).array() + 0.5) *
                                              (points.row(0).array() + points.row(1).array() + 0.5) *
                                              (points.row(0).array() + points.row(1).array() + 0.5) +
                                          time_value;

        // const Eigen::ArrayXd polynomial = (points.row(0).array() + points.row(1).array() + 0.5) + time_value;

        return polynomial;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int dimension,
                                                             const unsigned int id_domain,
                                                             const Eigen::MatrixXd &points,
                                                             const double &time_value) const
    {
        Eigen::VectorXd derivatives = (points.row(0).array() + points.row(1).array() + 0.5) *
                                      (points.row(0).array() + points.row(1).array() + 0.5) * 3.0;

        // Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), 1.0);

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }

    Eigen::VectorXd alpha(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    }

    Eigen::VectorXd beta(const unsigned int dimension, const unsigned int id_domain, const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    }
};
// ***************************************************************************
// struct Elliptic_Problem final : public I_Test
// {
//     Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
//     {
//         Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

//         domain.area = 1.0;

//         domain.vertices = Eigen::MatrixXd::Zero(3, 4);
//         domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
//         domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

//         domain.shape_type =
//         Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

//         return domain;
//     }

//     std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
//     {
//         return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
//                 {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
//                 {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
//     }

//     Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
//     {
//         const double k = 2.0;
//         return Eigen::VectorXd::Constant(points.cols(), k);
//     };

//     std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
//     {
//         return {Eigen::VectorXd::Constant(points.cols(), 0.0),
//                 Eigen::VectorXd::Constant(points.cols(), 0.0),
//                 Eigen::VectorXd::Constant(points.cols(), 0.0)};
//     };

//     Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
//     {
//         return 16.0 * std::numbers::pi * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
//                sin(2.0 * std::numbers::pi * points.row(1).array());
//     };

//     Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
//     {
//         if (marker != 1)
//             throw std::runtime_error("Unknown marker");

//         return exact_solution(points);
//     };

//     Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
//     {
//         switch (marker)
//         {
//         default:
//             throw std::runtime_error("Unknown marker");
//         }
//     }

//     Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
//     {
//         return sin(2.0 * std::numbers::pi * points.row(0).array()) * sin(2.0 * std::numbers::pi *
//         points.row(1).array());
//     };

//     std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
//     {
//         return {2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * points.row(0).array()) *
//                     sin(2.0 * std::numbers::pi * points.row(1).array()),
//                 2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
//                     cos(2.0 * std::numbers::pi * points.row(1).array()),
//                 Eigen::VectorXd::Zero(points.cols())};
//     }
// };
// ***************************************************************************
} // namespace test
} // namespace Parabolic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim

#endif
