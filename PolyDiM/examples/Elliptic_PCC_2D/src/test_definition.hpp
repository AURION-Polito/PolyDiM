#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Elliptic_Polynomial_Problem = 2, /// Test 1: S. Berrone, G. Teora, F. Vicini, "Improving high-order VEM stability on
                                     /// badly-shaped elements", doi: https://doi.org/10.1016/j.matcom.2023.10.003.
    SUPG_AdvDiff_Problem = 3 /// Test 1: M. Benedetto, S. Berrone, A. Borio, S. Pieraccini, S. Scialò, "Order preserving
                             /// SUPG stabilization for the Virtual Element formulation of advection-diffusion
                             /// problems", doi: https://doi.org/10.1016/j.cma.2016.07.043.
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const = 0;
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
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Elliptic_Polynomial_Problem final : public I_Test
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

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 4: // co-normal derivatives on the left
            return -16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct SUPG_AdvDiff_Problem final : public I_Test
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
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        const double k = 1.0e-06;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd &points) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 0.5),
                Eigen::VectorXd::Constant(points.cols(), -1.0 / 3.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        return -65536.0 / 729.0 * 1.0e-06 *
                   ((6.0 - 12.0 * points.row(0).array()) * points.row(0).array() * points.row(1).array() *
                        points.row(1).array() * points.row(1).array() * (1.0 - points.row(1).array()) +
                    (6.0 - 12.0 * points.row(1).array()) * points.row(1).array() * points.row(0).array() *
                        points.row(0).array() * points.row(0).array() * (1.0 - points.row(0).array())) +
               0.5 * 65536.0 / 729.0 * (3.0 - 4.0 * points.row(0).array()) * points.row(0).array() * points.row(0).array() *
                   points.row(1).array() * points.row(1).array() * points.row(1).array() * (1.0 - points.row(1).array()) -
               1.0 / 3.0 * 65536.0 / 729.0 * (3.0 - 4.0 * points.row(1).array()) * points.row(1).array() * points.row(1).array() *
                   points.row(0).array() * points.row(0).array() * points.row(0).array() * (1.0 - points.row(0).array());
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 65536.0 / 729.0 * points.row(1).array() * points.row(1).array() * points.row(1).array() *
               (1.0 - points.row(1).array()) * points.row(0).array() * points.row(0).array() * points.row(0).array() *
               (1.0 - points.row(0).array());
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return 65536.0 / 729.0 * points.row(1).array() * points.row(1).array() * points.row(1).array() *
               (1.0 - points.row(1).array()) * points.row(0).array() * points.row(0).array() * points.row(0).array() *
               (1.0 - points.row(0).array());
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {65536.0 / 729.0 * (3.0 - 4.0 * points.row(0).array()) * points.row(0).array() * points.row(0).array() *
                    points.row(1).array() * points.row(1).array() * points.row(1).array() * (1.0 - points.row(1).array()),
                65536.0 / 729.0 * (3.0 - 4.0 * points.row(1).array()) * points.row(1).array() * points.row(1).array() *
                    points.row(0).array() * points.row(0).array() * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
