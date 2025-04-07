#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Poisson_Polynomial_Problem = 2 /// Test 3: S. Berrone, G. Teora, F. Vicini, "Improving high-order VEM stability on
    /// badly-shaped elements", doi: https://doi.org/10.1016/j.matcom.2023.10.003.
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const = 0;
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
                {11, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {12, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {13, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {14, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {15, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {16, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {17, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {18, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {19, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {20, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {21, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {22, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {23, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {24, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {25, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {26, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}}};
    }

    Eigen::VectorXd diffusion_term(const Eigen::MatrixXd &points) const
    {
        return Eigen::VectorXd::Constant(points.cols(), 2.0);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 6.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), 2.0 * order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        switch (marker)
        {
        case 2: // co-normal derivatives on the left
            return -derivatives;
        case 4: // co-normal derivatives on the right
            return derivatives;
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, derivatives};
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
        const double k = 1.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
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
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 64.0 * points.row(2).array() * (1.0 - points.row(2).array()) * points.row(1).array() *
                   (1.0 - points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) +
               1.7;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        case 2: // bottom face
            return -64.0 * (1.0 - 2.0 * points.row(2).array()) * points.row(0).array() * (1.0 - points.row(0).array()) *
                   points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return 64.0 * points.row(2).array() * (1.0 - points.row(2).array()) * points.row(1).array() *
                   (1.0 - points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) +
               1.7;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {64.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()) *
                    points.row(2).array() * (1.0 - points.row(2).array()),
                64.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()) *
                    points.row(2).array() * (1.0 - points.row(2).array()),
                64.0 * (1.0 - 2.0 * points.row(2).array()) * points.row(0).array() * (1.0 - points.row(0).array()) *
                    points.row(1).array() * (1.0 - points.row(1).array())};
    }
};
// ***************************************************************************
} // namespace test
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim

#endif
