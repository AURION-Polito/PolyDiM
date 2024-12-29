#include <gtest/gtest.h>

#include "test_DOFsManager.hpp"

#include "test_FEM_Tetrahedron_PCC_3D_LocalSpace.hpp"
#include "test_FEM_Triangle_PCC_2D_LocalSpace.hpp"

#include "test_VEM_PCC_2D_Inertia_LocalSpace.hpp"
#include "test_VEM_PCC_2D_LocalSpace.hpp"
#include "test_VEM_PCC_2D_Ortho_LocalSpace.hpp"
#include "test_VEM_PCC_3D_Inertia_LocalSpace.hpp"
#include "test_VEM_PCC_3D_LocalSpace.hpp"
#include "test_VEM_PCC_3D_Ortho_LocalSpace.hpp"

#include "test_VEM_MCC_2D_LocalSpace.hpp"
#include "test_VEM_MCC_3D_LocalSpace.hpp"

#include "test_VEM_DF_PCC_2D_LocalSpace.hpp"
#include "test_VEM_DF_PCC_3D_LocalSpace.hpp"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
