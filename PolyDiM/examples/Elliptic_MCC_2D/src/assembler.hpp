#ifndef __assembler_H
#define __assembler_H

#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"
#include "DOFsManager.hpp"

#include "VEM_MCC_VelocityLocalSpace_Data.hpp"

namespace Elliptic_MCC_2D
{
template<typename VEM_LocalSpace_Type>
class Assembler final
{
public:

    struct Elliptic_MCC_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> neumannMatrixA;
        Gedim::Eigen_Array<> rightHandSide;
        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solutionNeumann;
    };

    struct VEM_Performance_Result final
    {
        struct Cell2D_Performance final
        {
            unsigned int NumBoundaryQuadraturePoints = 0;
            unsigned int NumInternalQuadraturePoints = 0;
            Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data Analysis;
        };

        std::vector<Cell2D_Performance> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        std::vector<Eigen::VectorXd> cell0Ds_numeric_pressure;
        std::vector<Eigen::VectorXd> cell0Ds_exact_pressure;

        Eigen::VectorXd cell2Ds_error_L2_pressure;
        Eigen::VectorXd cell2Ds_norm_L2_pressure;
        double error_L2_pressure;
        double norm_L2_pressure;
        Eigen::VectorXd cell2Ds_error_L2_velocity;
        Eigen::VectorXd cell2Ds_norm_L2_velocity;
        double error_L2_velocity;
        double norm_L2_velocity;

        double mesh_size;

        double residual_norm;
    };

private:

    void ComputeStrongTerm(const Gedim::GeometryUtilities& geometryUtilities,
                           const Gedim::MeshMatricesDAO& mesh,
                           const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                           const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                           const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                           const std::function<Eigen::VectorXd(const unsigned int,
                                                               const Eigen::MatrixXd&)>& strong_boundary_condition,
                           Elliptic_MCC_2D_Problem_Data& assembler_data) const;

    void ComputeWeakTerm(const unsigned int cell2DIndex,
                         const Gedim::MeshMatricesDAO& mesh,
                         const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry& polygon,
                         const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData &dofs_data,
                         const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                         const Polydim::VEM::MCC::VEM_MCC_VelocityLocalSpace_Data &local_space_data,
                         const std::function<Eigen::VectorXd(const unsigned int,
                                                             const Eigen::MatrixXd&)>& weak_boundary_condition,
                         Elliptic_MCC_2D_Problem_Data& assembler_data) const;

public:
    Elliptic_MCC_2D_Problem_Data Assemble(const Gedim::GeometryUtilities& geometryUtilities,
                                          const Gedim::MeshMatricesDAO& mesh,
                                          const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo>& mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData>& dofs_data,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                          const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& advection_term,
                                          const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& reaction_term,
                                          const std::function<std::array<Eigen::VectorXd, 9>(const Eigen::MatrixXd&)>& diffusion_term,
                                          const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& source_term,
                                          const std::function<Eigen::VectorXd(const unsigned int,
                                                                              const Eigen::MatrixXd&)>& strong_boundary_condition,
                                          const std::function<Eigen::VectorXd(const unsigned int,
                                                                              const Eigen::MatrixXd&)>& weak_boundary_condition) const;

    VEM_Performance_Result ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                 const Gedim::MeshMatricesDAO& mesh,
                                                 const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                 const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Gedim::GeometryUtilities& geometryUtilities,
                                         const Gedim::MeshMatricesDAO& mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                         const std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData>& dofs_data,
                                         const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                         const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                         const Elliptic_MCC_2D_Problem_Data& assembler_data,
                                         const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& exact_velocity,
                                         const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& exact_pressure) const;
};

}

#endif
