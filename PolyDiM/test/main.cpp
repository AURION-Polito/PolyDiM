// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include <gtest/gtest.h>

#include "test_DOFsManager.hpp"

#include "test_FEM_Hexahedron_PCC_3D_LocalSpace.hpp"
#include "test_FEM_PCC_1D_LocalSpace.hpp"
#include "test_FEM_Quadrilateral_PCC_2D_LocalSpace.hpp"
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
