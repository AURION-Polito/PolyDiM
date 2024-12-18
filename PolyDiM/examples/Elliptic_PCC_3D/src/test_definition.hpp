#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

#include <unordered_map>
#include <typeindex>

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
namespace test
{
// ***************************************************************************
struct Patch_Test final
{
    static unsigned int order;

    static Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain()
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0)<< 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1)<< 0.0, 0.0, 1.0, 1.0;

        return domain;
    }

    static std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info()
    {
        return {
            { 0, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } }
        };
    }

    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
    { return Eigen::VectorXd::Constant(points.cols(), 1.0); };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(),
                                                                2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return - source_term;
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    { throw std::runtime_error("Not supported"); }

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return
            {
                derivatives,
                derivatives,
                Eigen::VectorXd::Zero(points.cols())
            };
    }
};
unsigned int Patch_Test::order;
// ***************************************************************************
struct Poisson_Polynomial_Problem final
{
    static Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain()
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0)<< 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1)<< 0.0, 0.0, 1.0, 1.0;

        return domain;
    }

    static std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info()
    {
        return {
            { 0, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 }  },
            { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } },
            { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1 } },
            { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4 } }
        };
    }

    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
    {
        const double k = 1.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                       points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    {
        switch(marker)
        {
        case 2: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 4: // co-normal derivatives on the left
            return - 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                       points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    };

    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
    {
        return
            {
                16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())
            };
    }
};
// ***************************************************************************
// Centralized mapping function for IDs
unsigned int TestType(const std::type_index& type)
{
    static const std::unordered_map<std::type_index, unsigned int> typeToID = {
        {typeid(Polydim::examples::Elliptic_PCC_2D::test::Patch_Test), 1},
        {typeid(Polydim::examples::Elliptic_PCC_2D::test::Poisson_Polynomial_Problem), 2}
    };

    auto it = typeToID.find(type);
    if (it != typeToID.end()) {
        return it->second;
    } else {
        throw std::runtime_error("Class type not recognized.");
    }
}
// Helper template to get the ID for a specific class
template <typename T>
unsigned int TestType() {
    return TestType(typeid(T));
}
// ***************************************************************************
}
}
}
}

#endif
