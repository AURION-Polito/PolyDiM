#ifndef __test_definition_H
#define __test_definition_H

#include "CommonUtilities.hpp"
#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
namespace test
{
// ***************************************************************************
enum struct Test_Types
{
    Patch_Test = 1,
    StokesSinSin = 2,          /// Test 6.1 - Da Veiga Lovadina 2017
    Stokes_ZeroVelocity_1 = 3, /// Test 5.1 (a) - Da Veiga Lovadina 2018
    Stokes_ZeroVelocity_2 = 4, /// Test 5.1 (b) - Da Veiga Lovadina 2018
    Darcy = 5, /// Test 6.1 - Vacca 2018 "An H 1-conforming Virtual Element Methods for Darcy equations and Brinkman
    /// equations"
    Brinkman = 6, /// Test 6.2 - Vacca 2018 "An H 1-conforming Virtual Element Methods for Darcy equations and Brinkman
    /// equations"
    DarcyStokes_1 = 7,
    DarcyStokes_2 = 8 /// Test 5.6 - Da Veiga, Lovadina, Vacca, "VIRTUAL ELEMENTS FOR THE NAVIER--STOKES PROBLEM ON
    /// POLYGONAL MESHES", 2018, doi: 10.1137/17M1132811
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const = 0;
    virtual Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const = 0;
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

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

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 2.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 2; i++)
            result = result * polynomial;

        return {-4.0 * order * (order - 1) * result + (order - 1) * result + 0.8 * result * polynomial * polynomial,
                -4.0 * order * (order - 1) * result + (order - 1) * result + 0.6 * result * polynomial * polynomial,
                Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
            result = result * polynomial;

        return 2.0 * order * result;
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
        double mean = 0.0;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
        {
            result = result * polynomial;
            mean += Gedim::Utilities::BinomialCoefficient(order - 1.0, i) / ((i + 1.0) * (order - i));
        }

        mean += 1.0 / order;

        switch (marker)
        {
        case 2: // co-normal derivatives on the bottom
            return {2.0 * order * result, 2.0 * order * result - (result - mean), Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();
        double mean = 0.0;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 1; i++)
        {
            result = result * polynomial;
            mean += Gedim::Utilities::BinomialCoefficient(order - 1.0, i) / ((i + 1.0) * (order - i));
        }

        mean += 1.0 / order;

        result -= mean;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; i++)
            result = result * polynomial;

        return {result, result, Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
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
struct StokesSinSin final : public I_Test
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        std::vector<Eigen::VectorXd> laplacian(2);
        laplacian[0] = cos(points.row(1).array()) * sin(points.row(1).array()) *
                       (-sin(points.row(0).array()) * sin(points.row(0).array()) +
                        3.0 * cos(points.row(0).array()) * cos(points.row(0).array()));
        laplacian[1] = -cos(points.row(0).array()) * sin(points.row(0).array()) *
                       (-sin(points.row(1).array()) * sin(points.row(1).array()) +
                        3.0 * cos(points.row(1).array()) * cos(points.row(1).array()));

        std::vector<Eigen::VectorXd> pressure_derivatives(2);
        pressure_derivatives[0] = cos(points.row(0).array());
        pressure_derivatives[1] = -cos(points.row(1).array());

        return {-laplacian[0] + pressure_derivatives[0],
                -laplacian[1] + pressure_derivatives[1],
                Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::ArrayXd::Constant(points.cols(), 0.0);
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return {-0.5 * cos(points.row(0).array()) * cos(points.row(0).array()) * cos(points.row(1).array()) *
                    sin(points.row(1).array()),
                0.5 * cos(points.row(1).array()) * cos(points.row(1).array()) * cos(points.row(0).array()) *
                    sin(points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2:
            return {-0.5 * cos(points.row(0).array()) * cos(points.row(0).array()) *
                        (-sin(points.row(1).array()) * sin(points.row(1).array()) +
                         cos(points.row(1).array()) * cos(points.row(1).array())),
                    -cos(points.row(1).array()) * sin(points.row(1).array()) * cos(points.row(0).array()) *
                            sin(points.row(0).array()) -
                        sin(points.row(0).array()) - sin(points.row(1).array()),
                    Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return sin(points.row(0).array()) - sin(points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {-0.5 * cos(points.row(0).array()) * cos(points.row(0).array()) * cos(points.row(1).array()) *
                    sin(points.row(1).array()),
                0.5 * cos(points.row(1).array()) * cos(points.row(1).array()) * cos(points.row(0).array()) *
                    sin(points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {cos(points.row(0).array()) * sin(points.row(0).array()) * cos(points.row(1).array()) * sin(points.row(1).array()),
                -0.5 * cos(points.row(0).array()) * cos(points.row(0).array()) *
                    (-sin(points.row(1).array()) * sin(points.row(1).array()) +
                     cos(points.row(1).array()) * cos(points.row(1).array())),
                Eigen::VectorXd::Zero(points.cols()),
                0.5 * cos(points.row(1).array()) * cos(points.row(1).array()) *
                    (-sin(points.row(0).array()) * sin(points.row(0).array()) +
                     cos(points.row(0).array()) * cos(points.row(0).array())),
                -cos(points.row(1).array()) * sin(points.row(1).array()) * cos(points.row(0).array()) *
                    sin(points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Stokes_ZeroVelocity_1 final : public I_Test
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        return {3.0 * points.row(0).array() * points.row(0).array(),
                -3.0 * points.row(1).array() * points.row(1).array(),
                Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::ArrayXd::Constant(points.cols(), 0.0);
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
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(0).array() * points.row(0).array() -
               points.row(1).array() * points.row(1).array() * points.row(1).array();
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Stokes_ZeroVelocity_2 final : public I_Test
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        return {2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * points.row(0).array()) *
                    sin(2.0 * std::numbers::pi * points.row(1).array()),
                2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
                    cos(2.0 * std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::ArrayXd::Constant(points.cols(), 0.0);
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
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return sin(2.0 * std::numbers::pi * points.row(0).array()) * sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Darcy final : public I_Test
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 6}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 8}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 0.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return 2.0 * std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) *
               cos(std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Unknown marker");
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        const Eigen::VectorXd pressure =
            cos(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array());
        switch (marker)
        {
        case 2:
            return {Eigen::VectorXd::Zero(points.cols()), -pressure, Eigen::VectorXd::Zero(points.cols())};
        case 4:
            return {pressure, Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
        case 6:
            return {Eigen::VectorXd::Zero(points.cols()), pressure, Eigen::VectorXd::Zero(points.cols())};
        case 8:
            return {-pressure, Eigen::VectorXd::Zero(points.cols()), Eigen::VectorXd::Zero(points.cols())};
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return cos(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) *
                    cos(std::numbers::pi * points.row(1).array()),
                -std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) *
                    sin(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                -std::numbers::pi * std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) *
                    sin(std::numbers::pi * points.row(1).array()),
                std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) *
                    cos(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Brinkman final : public I_Test
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

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result = {
                                                                                                           {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                           {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                           {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result, result, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        return {(2.0 * std::numbers::pi * std::numbers::pi + 1.0) * sin(std::numbers::pi * points.row(0).array()) *
                        cos(std::numbers::pi * points.row(1).array()) +
                    2.0 * points.row(0).array() * points.row(1).array() * points.row(1).array(),
                (-2.0 * std::numbers::pi * std::numbers::pi - 1.0) * cos(std::numbers::pi * points.row(0).array()) *
                        sin(std::numbers::pi * points.row(1).array()) +
                    2.0 * points.row(1).array() * points.row(0).array() * points.row(0).array(),
                Eigen::VectorXd::Zero(points.cols())};
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Zero(points.cols());
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return {sin(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                -cos(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(0).array() * points.row(1).array() * points.row(1).array() - 1.0 / 9.0;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {sin(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                -cos(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {
                std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                -std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                std::numbers::pi * sin(std::numbers::pi * points.row(0).array()) * sin(std::numbers::pi * points.row(1).array()),
                -std::numbers::pi * cos(std::numbers::pi * points.row(0).array()) * cos(std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct DarcyStokes_1 final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 4.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 2.0, 2.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 2.0, 2.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return domain;
    }

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result1 = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                                                                                                            {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 12}},
                                                                                                            {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                                                                                                            {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 6}},
                                                                                                            {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 8}},
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 10}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result2 = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 14}},
                                                                                                            {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 16}},
                                                                                                            {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 18}},
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
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
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result1, result2, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        if (points(0, 0) < 1.0 - 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else if (points(1, 0) > 1.0 + 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 10.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 10.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else if (points(1, 0) < 1.0 - 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 2.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 2.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else
            throw std::runtime_error("not valid configuration");
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        if (points(0, 0) < 1.0 - 1.0e-12)
            return Eigen::VectorXd::Constant(points.cols(), 2.0);
        else if (points(0, 0) > 1.0 + 1.0e-12)
            return Eigen::VectorXd::Constant(points.cols(), 0.0);
        else
            throw std::runtime_error("not valid configuration");
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0).array();
        const Eigen::ArrayXd y = points.row(1).array();
        const Eigen::ArrayXd polynomial = x + y;

        if (points(0, 0) < 1.0 - 1.0e-12)
            return {Eigen::VectorXd::Ones(points.cols()),
                    Eigen::VectorXd::Ones(points.cols()),
                    Eigen::VectorXd::Zero(points.cols())};
        else if (points(1, 0) > 1.0 + 1.0e-12)
            return {Eigen::ArrayXd::Ones(points.cols()) + 10.0 * (1.0 - x),
                    Eigen::ArrayXd::Ones(points.cols()) + 10.0 * y,
                    Eigen::VectorXd::Zero(points.cols())};
        else if (points(1, 0) < 1.0 - 1.0e-12)
            return {Eigen::ArrayXd::Ones(points.cols()) + 2.0 * (1.0 - x),
                    Eigen::ArrayXd::Ones(points.cols()) + 2.0 * y,
                    Eigen::VectorXd::Zero(points.cols())};
        else
            throw std::runtime_error("not valid configuration");
    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Zero(points.cols());
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {

        switch (marker)
        {
        case 1:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = exact_velocity(points)[0];
            return result;
        }
        case 3:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = exact_velocity(points)[1];
            return result;
        }
        default:
            throw std::runtime_error("Unknown marker");
        }

    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0).array();
        const Eigen::ArrayXd y = points.row(1).array();
        const Eigen::ArrayXd polynomial = x + y;

        switch (marker)
        {
        case 2:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 4:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 6:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = polynomial;
            return result;
        }
        case 8:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = polynomial;
            return result;
        }
        case 10:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 12:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 14:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 16:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 18:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0).array();
        const Eigen::ArrayXd y = points.row(1).array();
        const Eigen::ArrayXd polynomial = x + y;

        bool left = false;
        for(unsigned int p = 0; p < points.cols(); p++)
        {
            if (points(0, p) < 1.0 - 1.0e-12)
            {
                left = true;
                break;
            }
        }

        if (left)
            return polynomial - 2.0;
        else
            return polynomial;
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0).array();
        const Eigen::ArrayXd y = points.row(1).array();

        return {1.0 - x, y, Eigen::VectorXd::Zero(points.cols())};

    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {-Eigen::VectorXd::Ones(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Ones(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct DarcyStokes_2 final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 4.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 2.0, 2.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 2.0, 2.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return domain;
    }

    std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 4> boundary_info() const
    {
        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result1 = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                                                                                                            {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 12}},
                                                                                                            {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                                                                                                            {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                                                                                                            {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 6}},
                                                                                                            {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 8}},
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 10}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> result2 = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {9, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {10, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 14}},
                                                                                                            {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
                                                                                                            {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 16}},
                                                                                                            {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                                                                                                            {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 18}},
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}}};

        std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> resultI = {
                                                                                                            {0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
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
                                                                                                            {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}}};

        return {result1, result2, resultI, resultI};
    }

    std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points) const
    {
        if (points(0, 0) < 1.0 - 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else if (points(1, 0) > 1.0 + 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 10.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 10.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else if (points(1, 0) < 1.0 - 1.0e-12)
            return {Eigen::VectorXd::Constant(points.cols(), 2.0),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0),
                    Eigen::VectorXd::Constant(points.cols(), 2.0),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Zero(points.cols()),
                    Eigen::VectorXd::Constant(points.cols(), 0.0)};
        else
            throw std::runtime_error("not valid configuration");
    };

    Eigen::VectorXd fluid_viscosity(const Eigen::MatrixXd &points) const
    {
        if (points(0, 0) < 1.0 - 1.0e-12)
            return Eigen::VectorXd::Constant(points.cols(), 1.0);
        else if (points(0, 0) > 1.0 + 1.0e-12)
            return Eigen::VectorXd::Constant(points.cols(), 0.0);
        else
            throw std::runtime_error("not valid configuration");
    };

    std::array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd &points) const
    {

        return {Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};

    };

    Eigen::VectorXd divergence_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Zero(points.cols());
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {

        switch (marker)
        {
        case 1:
        {
            std::array<Eigen::VectorXd, 3> result;

            if (points(1, 0) > 1.0 + 1.0e-12)
            {
                const Eigen::ArrayXd y = points.row(1);
                result[0] = -10.0 * (1.0 - y) * (2.0 - y);
            }
            else
                result[0] = Eigen::VectorXd::Zero(points.cols());

            return result;
        }
        case 3:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        default:
            throw std::runtime_error("Unknown marker");
        }

    }

    std::array<Eigen::VectorXd, 3> weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 4:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 6:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 8:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 10:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 12:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[0] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 14:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 16:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        case 18:
        {
            std::array<Eigen::VectorXd, 3> result;
            result[1] = Eigen::VectorXd::Zero(points.cols());
            return result;
        }
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not implemented method");
    }
};
// ***************************************************************************
} // namespace test
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
