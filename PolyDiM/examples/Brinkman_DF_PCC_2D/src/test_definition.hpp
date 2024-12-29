#ifndef __test_definition_H
#define __test_definition_H

#include "CommonUtilities.hpp"
#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

#include <unordered_map>

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
    Darcy = 5,                 /// Test 6.1 - Vacca 2018
    Brinkman = 6               /// Test 6.2 - Vacca 2018
};
// ***************************************************************************
struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
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

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
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
        return {2.0 * M_PI * cos(2.0 * M_PI * points.row(0).array()) * sin(2.0 * M_PI * points.row(1).array()),
                2.0 * M_PI * sin(2.0 * M_PI * points.row(0).array()) * cos(2.0 * M_PI * points.row(1).array()),
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
        return sin(2.0 * M_PI * points.row(0).array()) * sin(2.0 * M_PI * points.row(1).array());
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
        return 2.0 * M_PI * M_PI * cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return {M_PI * sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                M_PI * cos(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
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
        return cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd &points) const
    {
        return {M_PI * sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                M_PI * cos(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {M_PI * M_PI * cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                -M_PI * M_PI * sin(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                -M_PI * M_PI * sin(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                M_PI * M_PI * cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
struct Brinkman final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

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
        return {(2.0 * M_PI * M_PI + 1.0) * sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()) +
                    2.0 * points.row(0).array() * points.row(1).array() * points.row(1).array(),
                (-2.0 * M_PI * M_PI - 1.0) * sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()) +
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

        return {sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                -cos(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
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
        return {sin(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                -cos(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }

    std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd &points) const
    {
        return {M_PI * cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                -M_PI * sin(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                M_PI * sin(M_PI + points.row(0).array()) * sin(M_PI + points.row(1).array()),
                -M_PI * cos(M_PI + points.row(0).array()) * cos(M_PI + points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};

// ***************************************************************************
} // namespace test
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
