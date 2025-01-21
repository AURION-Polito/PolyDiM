#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

#include <numbers>
#include <unordered_map>

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
namespace test
{
// ***************************************************************************
enum struct Test_Types
{
    Patch_Test = 1,
    Stokes = 2,             /// Example 1 - Da Veiga Dassi 2020
    Stokes_Benchmark_1 = 3, /// Example 3.a - Da Veiga Dassi 2020
    Stokes_Benchmark_2 = 4  /// Example 3.b - Da Veiga Dassi 2020
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker,
                                                                     const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker,
                                                                   const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const = 0;
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

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 2; i++)
            result = result * polynomial;

        return {-3.0 * order * (order - 1) * result - (order - 1) * result,
                -3.0 * order * (order - 1) * result - (order - 1) * result,
                6.0 * order * (order - 1) * result - (order - 1) * result};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return {result, result, -2.0 * result};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("not implemented method");
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array();
        double mean = 0.0;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
        {
            result = result * polynomial;

            for (int j = 0; j <= order - 1 - i; j++)
                mean += std::tgamma(order) / (std::tgamma(i + 1) * std::tgamma(j + 1) * std::tgamma(order - i - j) *
                                              (i + 1.0) * (j + 1.0) * (order - i - j));
        }

        mean += 1.0 / order;

        result -= mean;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return {result, result, -2.0 * result};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
            result = result * polynomial;

        result *= order;

        return {result, result, result, result, result, result, -2.0 * result, -2.0 * result, -2.0 * result};
    }
};
// ***************************************************************************
struct Stokes final : public I_Test
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

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        std::vector<Eigen::VectorXd> derivatesTermValues(3, Eigen::VectorXd::Zero(points.cols()));
        derivatesTermValues[0] = std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * x) *
                                 cos(std::numbers::pi * y) * cos(std::numbers::pi * z);
        derivatesTermValues[1] = std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * x) *
                                 sin(std::numbers::pi * y) * cos(std::numbers::pi * z);
        derivatesTermValues[2] = std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * x) *
                                 cos(std::numbers::pi * y) * sin(std::numbers::pi * z);

        std::vector<Eigen::VectorXd> laplacianVelocityTermValues(3, Eigen::VectorXd::Zero(points.cols()));
        laplacianVelocityTermValues[0] = -3.0 * std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * x) *
                                         cos(std::numbers::pi * y) * cos(std::numbers::pi * z);
        laplacianVelocityTermValues[1] = -3.0 * std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * x) *
                                         sin(std::numbers::pi * y) * cos(std::numbers::pi * z);
        laplacianVelocityTermValues[2] = 6.0 * std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * x) *
                                         cos(std::numbers::pi * y) * sin(std::numbers::pi * z);

        return {-laplacianVelocityTermValues[0] - derivatesTermValues[0],
                -laplacianVelocityTermValues[1] - derivatesTermValues[1],
                -laplacianVelocityTermValues[2] - derivatesTermValues[2]};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        return {sin(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z),
                cos(std::numbers::pi * x) * sin(std::numbers::pi * y) * cos(std::numbers::pi * z),
                -2.0 * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * sin(std::numbers::pi * z)};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("not implemented method");
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        return -std::numbers::pi * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z);
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        return {sin(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z),
                cos(std::numbers::pi * x) * sin(std::numbers::pi * y) * cos(std::numbers::pi * z),
                -2.0 * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * sin(std::numbers::pi * z)};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        return {std::numbers::pi * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z),
                -std::numbers::pi * sin(std::numbers::pi * x) * sin(std::numbers::pi * y) * cos(std::numbers::pi * z),
                -std::numbers::pi * sin(std::numbers::pi * x) * cos(std::numbers::pi * y) * sin(std::numbers::pi * z),
                -std::numbers::pi * sin(std::numbers::pi * x) * sin(std::numbers::pi * y) * cos(std::numbers::pi * z),
                std::numbers::pi * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z),
                -std::numbers::pi * cos(std::numbers::pi * x) * sin(std::numbers::pi * y) * sin(std::numbers::pi * z),
                2.0 * std::numbers::pi * sin(std::numbers::pi * x) * cos(std::numbers::pi * y) * sin(std::numbers::pi * z),
                2.0 * std::numbers::pi * cos(std::numbers::pi * x) * sin(std::numbers::pi * y) * sin(std::numbers::pi * z),
                -2.0 * std::numbers::pi * cos(std::numbers::pi * x) * cos(std::numbers::pi * y) * cos(std::numbers::pi * z)};
    }
};
// ***************************************************************************
struct Stokes_Benchmark_1 final : public I_Test
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

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        std::vector<Eigen::VectorXd> derivativesPressureValues(3, Eigen::VectorXd::Zero(points.cols()));
        std::vector<Eigen::VectorXd> laplacianVelocityTermValues(3, Eigen::VectorXd::Zero(points.cols()));

        if (order == 2)
        {
            laplacianVelocityTermValues[2] = Eigen::VectorXd::Constant(points.cols(), -2.0 * order * (order - 1.0));

            derivativesPressureValues[0] = order * x * y + z * z;
            derivativesPressureValues[1] = order * y * z + x * x;
            derivativesPressureValues[2] = order * z * x + y * y;
        }
        else
        {
            Eigen::ArrayXd zkm3 = Eigen::VectorXd::Ones(points.cols());
            Eigen::ArrayXd xkm3 = Eigen::VectorXd::Ones(points.cols());
            Eigen::ArrayXd ykm3 = Eigen::VectorXd::Ones(points.cols());
            for (unsigned int i = 0; i < order - 3; i++)
            {
                zkm3 = zkm3 * z;
                xkm3 = xkm3 * x;
                ykm3 = ykm3 * y;
            }

            const Eigen::ArrayXd zkm1 = zkm3 * z * z;
            const Eigen::ArrayXd xkm1 = xkm3 * x * x;
            const Eigen::ArrayXd ykm1 = ykm3 * y * y;

            derivativesPressureValues[0] = order * xkm1 * y + zkm1 * z;
            derivativesPressureValues[1] = order * ykm1 * z + xkm1 * x;
            derivativesPressureValues[2] = order * zkm1 * x + ykm1 * y;

            laplacianVelocityTermValues[0] = order * (order - 1.0) * (order - 2) * x * zkm3;
            laplacianVelocityTermValues[1] = order * (order - 1.0) * (order - 2) * y * zkm3;
            laplacianVelocityTermValues[2] =
                (2.0 - order) * order * (order - 1.0) * (x * xkm3 + y * ykm3) - 2.0 * order * (order - 1.0) * zkm3 * z;
        }

        return {-laplacianVelocityTermValues[0] - derivativesPressureValues[0],
                -laplacianVelocityTermValues[1] - derivativesPressureValues[1],
                -laplacianVelocityTermValues[2] - derivativesPressureValues[2]};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm1 = Eigen::VectorXd::Ones(points.cols());
        for (unsigned int i = 0; i < order - 1; i++)
        {
            zkm1 = zkm1 * z;
            xkm1 = xkm1 * x;
            ykm1 = ykm1 * y;
        }

        return {order * x * zkm1, order * y * zkm1, (2.0 - order) * xkm1 * x + (2.0 - order) * ykm1 * y - 2.0 * zkm1 * z};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("not implemented method");
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zk = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xk = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd yk = Eigen::VectorXd::Ones(points.cols());
        for (unsigned int i = 0; i < order; i++)
        {
            zk = zk * z;
            xk = xk * x;
            yk = yk * y;
        }

        return xk * y + yk * z + zk * x - 3.0 / (2.0 * (order + 1.0));
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm1 = Eigen::VectorXd::Ones(points.cols());
        for (unsigned int i = 0; i < order - 1; i++)
        {
            zkm1 = zkm1 * z;
            xkm1 = xkm1 * x;
            ykm1 = ykm1 * y;
        }

        return {order * x * zkm1, order * y * zkm1, (2.0 - order) * xkm1 * x + (2.0 - order) * ykm1 * y - 2.0 * zkm1 * z};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm2 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm2 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm2 = Eigen::VectorXd::Ones(points.cols());

        for (unsigned int i = 0; i < order - 2; i++)
        {
            zkm2 = zkm2 * z;
            xkm2 = xkm2 * x;
            ykm2 = ykm2 * y;
        }

        return {order * zkm2 * z,
                Eigen::VectorXd::Zero(points.cols()),
                order * (order - 1) * x * zkm2,
                Eigen::VectorXd::Zero(points.cols()),
                order * zkm2 * z,
                order * (order - 1) * y * zkm2,
                (2.0 - order) * order * xkm2 * x,
                (2.0 - order) * order * ykm2 * y,
                -2.0 * order * zkm2 * z};
    }
};
// ***************************************************************************
struct Stokes_Benchmark_2 final : public I_Test
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

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        std::vector<Eigen::VectorXd> derivativesPressureValues(3, Eigen::VectorXd::Zero(points.cols()));
        derivativesPressureValues[0] = 2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * x) *
                                       sin(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * z);
        derivativesPressureValues[1] = 2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * x) *
                                       cos(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * z);
        derivativesPressureValues[2] = 2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * x) *
                                       sin(2.0 * std::numbers::pi * y) * cos(2.0 * std::numbers::pi * z);

        std::vector<Eigen::VectorXd> laplacianVelocityTermValues(3, Eigen::VectorXd::Zero(points.cols()));

        if (order == 2)
        {
            laplacianVelocityTermValues[2] = Eigen::VectorXd::Constant(points.cols(), -2.0 * order * (order - 1.0));
        }
        else
        {
            Eigen::ArrayXd zkm3 = Eigen::VectorXd::Ones(points.cols());
            Eigen::ArrayXd xkm3 = Eigen::VectorXd::Ones(points.cols());
            Eigen::ArrayXd ykm3 = Eigen::VectorXd::Ones(points.cols());
            for (unsigned int i = 0; i < order - 3; i++)
            {
                zkm3 = zkm3 * z;
                xkm3 = xkm3 * x;
                ykm3 = ykm3 * y;
            }

            laplacianVelocityTermValues[0] = order * (order - 1.0) * (order - 2) * x * zkm3;
            laplacianVelocityTermValues[1] = order * (order - 1.0) * (order - 2) * y * zkm3;
            laplacianVelocityTermValues[2] =
                (2.0 - order) * order * (order - 1.0) * (x * xkm3 + y * ykm3) - 2.0 * order * (order - 1.0) * zkm3 * z;
        }

        return {-laplacianVelocityTermValues[0] - derivativesPressureValues[0],
                -laplacianVelocityTermValues[1] - derivativesPressureValues[1],
                -laplacianVelocityTermValues[2] - derivativesPressureValues[2]};
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm1 = Eigen::VectorXd::Ones(points.cols());
        for (unsigned int i = 0; i < order - 1; i++)
        {
            zkm1 = zkm1 * z;
            xkm1 = xkm1 * x;
            ykm1 = ykm1 * y;
        }

        return {order * x * zkm1, order * y * zkm1, (2.0 - order) * xkm1 * x + (2.0 - order) * ykm1 * y - 2.0 * zkm1 * z};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("not implemented method");
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        return sin(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * y) * sin(2.0 * std::numbers::pi * z);
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm1 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm1 = Eigen::VectorXd::Ones(points.cols());
        for (unsigned int i = 0; i < order - 1; i++)
        {
            zkm1 = zkm1 * z;
            xkm1 = xkm1 * x;
            ykm1 = ykm1 * y;
        }

        return {order * x * zkm1, order * y * zkm1, (2.0 - order) * xkm1 * x + (2.0 - order) * ykm1 * y - 2.0 * zkm1 * z};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);
        Eigen::ArrayXd zkm2 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd xkm2 = Eigen::VectorXd::Ones(points.cols());
        Eigen::ArrayXd ykm2 = Eigen::VectorXd::Ones(points.cols());

        for (unsigned int i = 0; i < order - 2; i++)
        {
            zkm2 = zkm2 * z;
            xkm2 = xkm2 * x;
            ykm2 = ykm2 * y;
        }

        return {order * zkm2 * z,
                Eigen::VectorXd::Zero(points.cols()),
                order * (order - 1) * x * zkm2,
                Eigen::VectorXd::Zero(points.cols()),
                order * zkm2 * z,
                order * (order - 1) * y * zkm2,
                (2.0 - order) * order * xkm2 * x,
                (2.0 - order) * order * ykm2 * y,
                -2.0 * order * zkm2 * z};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
