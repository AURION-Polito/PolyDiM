#ifndef __local_space_H
#define __local_space_H

#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_Creator.hpp"
#include "program_configuration.hpp"

namespace Polydim
{
  namespace examples
  {
    namespace Elliptic_PCC_2D
    {
      namespace local_space
      {
        struct ReferenceElement_Data final
        {
            Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement;
            Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data FEM_ReferenceElement;
        };

        struct LocalSpace_Data final
        {
            Program_configuration::MethodTypes Method_Type;

            VEM::PCC::VEM_PCC_2D_LocalSpace_Types VEM_Type;
            Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry VEM_Geometry;
            std::unique_ptr<VEM::PCC::I_VEM_PCC_2D_LocalSpace> VEM_LocalSpace;
            Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data VEM_LocalSpace_Data;

            Polydim::FEM::PCC::FEM_Triangle_PCC_2D_Polygon_Geometry FEM_Geometry;
            std::unique_ptr<FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace> FEM_LocalSpace;
            Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data FEM_LocalSpace_Data;
        };

        LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_PCC_2D::Program_configuration &config,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const unsigned int cell2D_index,
                                         const ReferenceElement_Data& reference_element_data);

        Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data& reference_element_data,
                                             const LocalSpace_Data& local_space_data);


        std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(const ReferenceElement_Data& reference_element_data,
                                                                    const LocalSpace_Data& local_space_data);

        Eigen::MatrixXd StabilizationMatrix(const LocalSpace_Data& local_space_data);


        Gedim::Quadrature::QuadratureData InternalQuadrature(const LocalSpace_Data& local_space_data);

        unsigned int Size(const LocalSpace_Data& local_space_data);

      };
    } // namespace Elliptic_PCC_2D
  } // namespace examples
} // namespace Polydim

#endif
