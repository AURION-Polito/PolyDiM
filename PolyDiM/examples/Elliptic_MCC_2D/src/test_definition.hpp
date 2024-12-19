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
namespace Elliptic_MCC_2D
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
            { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0 } },
            { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } },
            { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } },
            { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } },
            { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2 } }
        };
    }

    static std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Constant(points.cols(), -0.2),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static Eigen::VectorXd reaction_term(const Eigen::MatrixXd& points)
    {
        return points.row(0).array() * points.row(1).array();
    }

    static std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 2.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };

    static std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 0.6),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
        Eigen::ArrayXd second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if(order > 1)
        {
            second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for(int i = 0; i < order - 2; i++)
                second_derivatives = second_derivatives * polynomial;

            solution = second_derivatives * polynomial * polynomial;
            second_derivatives *= order * (order - 1);
        }
        else if(order == 1)
            solution = polynomial;

        return - 3.0 * second_derivatives + points.row(1).array().transpose() * points.row(0).array().transpose() * solution;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
        switch(marker)
        {
        case 1: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 3: // co-normal derivatives on the left
            return - 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    static Eigen::VectorXd exact_pressure(const Eigen::MatrixXd& points)
    {

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    static std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd& points)
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if(order > 0)
        {
            derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for(int i = 0; i < order - 1; i++)
                derivatives = derivatives * polynomial;

            solution = derivatives * polynomial;
            derivatives *= order;
        }

        return
            {
                -derivatives + solution,
                -2.0 * derivatives - solution,
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
            { 1, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 0 } },
            { 2, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 0 } },
            { 3, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 0 } },
            { 4, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 0 } },
            { 5, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 2 } },
            { 6, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 1 } },
            { 7, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 2 } },
            { 8, { Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 3 } }
        };
    }

    static std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static Eigen::VectorXd reaction_term(const Eigen::MatrixXd& points)
    {
        return Eigen::VectorXd::Zero(points.cols());
    }

    static std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };

    static std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };


    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
        switch(marker)
        {
        case 1: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 3: // co-normal derivatives on the left
            return - 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                       points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    }

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                       points.row(0).array() * (1.0 - points.row(0).array())) + 1.1;
    };

    static std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd& points)
    {
        return
            {
                -16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                -16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())
            };
    }
};
// ***************************************************************************
// Centralized mapping function for IDs
unsigned int TestType(const std::type_index& type)
{
    static const std::unordered_map<std::type_index, unsigned int> typeToID = {
        {typeid(Polydim::examples::Elliptic_MCC_2D::test::Patch_Test), 1},
        {typeid(Polydim::examples::Elliptic_MCC_2D::test::Poisson_Polynomial_Problem), 2}
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
