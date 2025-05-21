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
namespace Elastic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    LinearElasticity = 2,
    LinearElasticity_Beam = 3,
    LinearElasticity_CooksMembrane = 4 /// Test 4 in Artioli, E., Beirão da Veiga, L., Lovadina, C., Sacco, E.,
                                       /// "Arbitrary order 2D virtual elements for polygonal meshes: part I, elastic
                                       /// problem", doi:  - 10.1007/s00466-017-1404-5.
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
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}}};
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

        const double mu = 1.0;
        const double lambda = 1.0;

        const Eigen::VectorXd u_1_x = order * result;
        const Eigen::VectorXd u_1_y = order * result;

        const Eigen::VectorXd u_2_x = order * result;
        const Eigen::VectorXd u_2_y = order * result;
        switch (marker)
        {
        case 2: // co-normal derivatives on the bottom
            return {-mu * (u_1_y + u_2_x), -((lambda + 2.0 * mu) * u_2_y + lambda * u_1_x), Eigen::VectorXd::Zero(points.cols())};
        case 4: // co-normal derivatives on the left
            return {-((lambda + 2.0 * mu) * u_1_x + lambda * u_1_y), -mu * (u_1_y + u_2_x), Eigen::VectorXd::Zero(points.cols())};
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
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
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
        const double mu = 1.5;
        const double lambda = 3.0;

        const Eigen::VectorXd u_1_x =
            16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        const Eigen::VectorXd u_1_y =
            16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array());

        const Eigen::VectorXd u_2_x =
            5.0 * 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        const Eigen::VectorXd u_2_y =
            5.0 * 16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array());
        switch (marker)
        {
        case 2: // co-normal derivatives on the top
            return {mu * (u_1_y + u_2_x), ((lambda + 2.0 * mu) * u_2_y + lambda * u_1_x), Eigen::VectorXd::Zero(points.cols())};
        case 4: // co-normal derivatives on the right
            return {((lambda + 2.0 * mu) * u_1_x + lambda * u_1_y), mu * (u_1_y + u_2_x), Eigen::VectorXd::Zero(points.cols())};
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
struct LinearElasticity_Beam final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 25.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 25.0, 25.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 6}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {

        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols())};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2: // co-normal derivatives on the bottom
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        case 4: // co-normal derivatives on the right
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        case 6: // co-normal derivatives on the top
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    std::array<Eigen::VectorXd, 2> lame_coefficients(const Eigen::MatrixXd &points) const
    {
        const double E = 1.0e05;
        const double nu = 0.3;
        const double mu = E * (2.0 * (1.0 + nu));
        const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        return {Eigen::ArrayXd::Constant(points.cols(), mu), Eigen::ArrayXd::Constant(points.cols(), lambda)};
    }

    std::array<Eigen::VectorXd, 3> exact_displacement(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_displacement(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }
};
// ***************************************************************************
struct LinearElasticity_CooksMembrane final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1408.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 48.0, 48.0, 0.0;
        domain.vertices.row(1) << 0.0, 44.0, 60.0, 44.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 6}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {

        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols())};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2: // co-normal derivatives on the bottom
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        case 4: // co-normal derivatives on the right
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 6.25),
                    Eigen::VectorXd::Zero(points.cols())};
        case 6: // co-normal derivatives on the top
            return {Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    std::array<Eigen::VectorXd, 2> lame_coefficients(const Eigen::MatrixXd &points) const
    {
        const double E = 70;
        const double nu = 1.0 / 3.0;
        const double mu = E * (2.0 * (1.0 + nu));
        const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        return {Eigen::ArrayXd::Constant(points.cols(), mu), Eigen::ArrayXd::Constant(points.cols(), lambda)};
    }

    std::array<Eigen::VectorXd, 3> exact_displacement(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_displacement(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
