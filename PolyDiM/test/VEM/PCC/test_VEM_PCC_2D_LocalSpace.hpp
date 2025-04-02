#ifndef __TEST_VEM_PCC_2D_LocalSpace_H
#define __TEST_VEM_PCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numbers>

#include "GeometryUtilities.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{

void Test_VEM_PCC_2D_Export_Dofs(const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry &polygon,
                                 const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &vem_reference_element,
                                 const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polygon.Vertices.cols();
    const unsigned int num_dofs = num_vertices * (vem_reference_element.Order + 1) + vem_reference_element.NumDofs2D;

    Eigen::MatrixXd dofs_coordinate = Eigen::MatrixXd::Zero(3, num_dofs);
    std::vector<double> dof_global_index_values(num_dofs);
    std::vector<double> dof_cell_index_values(num_dofs);
    std::vector<double> dof_dimension_values(num_dofs);

    unsigned int id_dofs = 0;
    if (vem_reference_element.NumDofs0D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs0D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 0;
                dofs_coordinate.col(id_dofs) = polygon.Vertices.col(c);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs1D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            const std::vector<double> local_edge_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs1D, 0.0, 1.0, false);
            const Eigen::Vector3d edge_origin = polygon.Vertices.col(c);
            const Eigen::Vector3d edge_tangent = polygon.Vertices.col((c + 1) % num_vertices) - edge_origin;

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs1D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 1;
                dofs_coordinate.col(id_dofs) = edge_origin + local_edge_coordinates[loc_i] * edge_tangent;
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs2D != 0)
    {

        const auto local_polygon_coordinates =
            geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = polygon.Centroid;
        const auto polygonCentroidEdgesDistance =
            geometryUtilities.PolygonCentroidEdgesDistance(polygon.Vertices, polygon.Centroid, polygon.EdgesNormal);

        double circle_diameter = 0.0;
        if (vem_reference_element.NumDofs2D > 1)
            circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs2D; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 2;
            dofs_coordinate.col(id_dofs) =
                polygon_centroid +
                circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  0.0);
            dof_global_index_values[id_dofs] = id_dofs;
            id_dofs++;
        }
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoints(dofs_coordinate,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values.size()),
                             dof_dimension_values.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values.size()),
                             dof_cell_index_values.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_global_index_values.size()),
                             dof_global_index_values.data()}});

        exporter.Export(exportVtuFolder + "/dofs.vtu");
    }
}

struct Test_VEM_PCC_2D_Polygon_Geometry final
{
    Eigen::MatrixXd Vertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::Matrix3d> TriangulationVertices;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;
};

Test_VEM_PCC_2D_Polygon_Geometry Test_VEM_PCC_2D_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_PCC_2D_Polygon_Geometry result;

    Eigen::MatrixXd &polygonVertices = result.Vertices;
    polygonVertices.setZero(3, 8);
    polygonVertices.row(0) << 2.8000000000000004e-03, 3.7200000000000004e-02, 3.6900000000000002e-02, 3.3300000000000003e-02,
        2.5600000000000001e-02, 1.4400000000000000e-02, 6.7000000000000002e-03, 3.0999999999999999e-03;
    polygonVertices.row(1) << 2.1499999999999998e-02, 2.1499999999999998e-02, 2.9500000000000002e-02, 4.0599999999999997e-02,
        5.0700000000000002e-02, 5.0700000000000002e-02, 4.0599999999999997e-02, 2.9500000000000002e-02;

    result.Measure = 7.9890999999999990e-04;
    result.Diameter = 3.7046997179258676e-02;
    result.Centroid = Eigen::Vector3d(2.0000000000000004e-02, 3.4061346918509809e-02, 0.0);
    result.EdgesLength.resize(8);
    result.EdgesLength << 3.4400000000000000e-02, 8.0056230238501787e-03, 1.1669190203266030e-02, 1.2700393694685222e-02,
        1.1200000000000002e-02, 1.2700393694685220e-02, 1.1669190203266030e-02, 8.0056230238501769e-03;

    result.EdgesDirection.resize(8, true);

    Eigen::MatrixXd &edgeTangents = result.EdgesTangent;
    edgeTangents.setZero(3, 8);
    edgeTangents.row(0) << 3.4400000000000000e-02, -3.0000000000000165e-04, -3.5999999999999990e-03, -7.7000000000000020e-03,
        -1.1200000000000002e-02, -7.6999999999999994e-03, -3.6000000000000003e-03, -2.9999999999999949e-04;
    edgeTangents.row(1) << 0.0000000000000000e+00, 8.0000000000000036e-03, 1.1099999999999995e-02, 1.0100000000000005e-02,
        0.0000000000000000e+00, -1.0100000000000005e-02, -1.1099999999999995e-02, -8.0000000000000036e-03;

    Eigen::MatrixXd &edgeNormals = result.EdgesNormal;
    edgeNormals.setZero(3, 8);
    edgeNormals.row(0) << 0.0000000000000000e+00, 9.9929761570918052e-01, 9.5122281894876248e-01, 7.9525093810490199e-01,
        0.0000000000000000e+00, -7.9525093810490211e-01, -9.5122281894876248e-01, -9.9929761570918074e-01;
    edgeNormals.row(1) << -1.0000000000000000e+00, 3.7473660589094460e-02, 3.0850469803743652e-01, 6.0628041815918254e-01,
        1.0000000000000000e+00, 6.0628041815918243e-01, 3.0850469803743663e-01, 3.7473660589094196e-02;

    // Map triangle points
    std::vector<unsigned int> polygonTriangulation = geometry_utilities.PolygonTriangulationByFirstVertex(polygonVertices);
    result.TriangulationVertices = geometry_utilities.ExtractTriangulationPoints(polygonVertices, polygonTriangulation);

    return result;
}

std::vector<Eigen::MatrixXd> Test_VEM_PCC_2D_RefPiNabla()
{
    std::vector<Eigen::MatrixXd> refPiNabla(3);
    refPiNabla[0] = (Eigen::MatrixXd(3, 8) << 1.7949235576524070e-01,
                     1.7949235576524064e-01,
                     9.0593575113167221e-02,
                     1.1461049108920020e-01,
                     1.1530357803239182e-01,
                     1.1530357803239186e-01,
                     1.1461049108920031e-01,
                     9.0593575113167318e-02,
                     -1.8548771290512855e-01,
                     1.8548771290512855e-01,
                     4.4285191456099415e-01,
                     4.9154243919859042e-01,
                     2.3417823754272479e-01,
                     -2.3417823754272479e-01,
                     -4.9154243919859042e-01,
                     -4.4285191456099415e-01,
                     -7.9064137625811004e-01,
                     -7.9064137625811004e-01,
                     9.0425260041250136e-02,
                     2.6200139447849402e-01,
                     4.3821472173836606e-01,
                     4.3821472173836606e-01,
                     2.6200139447849397e-01,
                     9.0425260041250122e-02)
                        .finished();
    refPiNabla[1] = (Eigen::MatrixXd(6, 17) << -5.9032716442944118e-02,
                     -5.9032716442944153e-02,
                     -3.1814258074315344e-02,
                     -3.7117946099565167e-02,
                     -3.8701746049843819e-02,
                     -3.8701746049843756e-02,
                     -3.7117946099565056e-02,
                     -3.1814258074315253e-02,
                     -1.8029161983066155e-01,
                     -5.5839245941115022e-02,
                     -7.1417786356146359e-02,
                     -7.7053998042114324e-02,
                     -7.7752986157260950e-02,
                     -7.7053998042114102e-02,
                     -7.1417786356146110e-02,
                     -5.5839245941114890e-02,
                     2.0000000000000200e+00,
                     -6.1829237635043124e-02,
                     6.1829237635042507e-02,
                     1.4761730485366395e-01,
                     1.6384747973286276e-01,
                     7.8059412514241422e-02,
                     -7.8059412514241686e-02,
                     -1.6384747973286398e-01,
                     -1.4761730485366537e-01,
                     5.5772846884634458e-32,
                     2.4731695054017006e-01,
                     3.4315226887448586e-01,
                     3.1223765005696535e-01,
                     1.0640911219645034e-32,
                     -3.1223765005696708e-01,
                     -3.4315226887448880e-01,
                     -2.4731695054017244e-01,
                     1.0660111602432053e-14,
                     -2.6354712541937070e-01,
                     -2.6354712541936903e-01,
                     3.0141753347083387e-02,
                     8.7333798159497600e-02,
                     1.4607157391278827e-01,
                     1.4607157391278894e-01,
                     8.7333798159498266e-02,
                     3.0141753347083331e-02,
                     -1.0634628873227359e+00,
                     9.2743856452567764e-03,
                     1.1129262774307669e-01,
                     2.3804256489491385e-01,
                     3.4624373075623965e-01,
                     2.3804256489491557e-01,
                     1.1129262774307733e-01,
                     9.2743856452560964e-03,
                     -8.5929004265210054e-17,
                     2.6718593835839033e-01,
                     2.6718593835838639e-01,
                     6.2678015292662825e-01,
                     5.4749816264194151e-01,
                     1.0982584791824485e-01,
                     1.0982584791824421e-01,
                     5.4749816264194151e-01,
                     6.2678015292662959e-01,
                     1.6991856616451224e-29,
                     1.0594233137233753e+00,
                     1.3018324165188786e+00,
                     7.4132447344815122e-01,
                     3.2418792981214800e-30,
                     7.4132447344815000e-01,
                     1.3018324165188797e+00,
                     1.0594233137233788e+00,
                     -9.3077406110712175e+00,
                     1.4516809119745731e+00,
                     -1.4516809119745637e+00,
                     -4.4821215412527704e-02,
                     6.1047707273402352e-01,
                     5.7874122836507202e-01,
                     -5.7874122836506936e-01,
                     -6.1047707273402030e-01,
                     4.4821215412530438e-02,
                     1.4420310806520018e-14,
                     -5.3566206910091196e-01,
                     5.5221370921984014e-01,
                     1.6043074655215361e+00,
                     2.7512536229191822e-15,
                     -1.6043074655215299e+00,
                     -5.5221370921983404e-01,
                     5.3566206910091752e-01,
                     -5.2468527856951727e-14,
                     9.9266583016712195e-01,
                     9.9266583016711785e-01,
                     -4.1225904694565073e-02,
                     1.7122975012678429e-01,
                     7.2877330250167782e-01,
                     7.2877330250167605e-01,
                     1.7122975012678252e-01,
                     -4.1225904694564934e-02,
                     4.0055958425511937e+00,
                     -2.3808707817444753e-02,
                     3.2992779160660325e-02,
                     8.2717264209975117e-01,
                     1.7274626429669380e+00,
                     8.2717264209974628e-01,
                     3.2992779160658771e-02,
                     -2.3808707817443070e-02,
                     -1.1108657868606089e+01)
                        .finished();
    refPiNabla[2] = (Eigen::MatrixXd(10, 27) << -2.6878028213698588e-02,
                     -2.6878028213698554e-02,
                     -1.5709627849840677e-02,
                     -2.2330918637009884e-02,
                     -1.8758817543076192e-02,
                     -1.8758817543076199e-02,
                     -2.2330918637009891e-02,
                     -1.5709627849840698e-02,
                     -1.2909059161907299e-01,
                     -1.2909059161907296e-01,
                     -2.9223295828573227e-02,
                     -3.1979490931262330e-02,
                     -4.7982595738891826e-02,
                     -5.1842645119807226e-02,
                     -5.7046739837807017e-02,
                     -5.3622550984169126e-02,
                     -4.2558373990761920e-02,
                     -4.2558373990761934e-02,
                     -5.3622550984169140e-02,
                     -5.7046739837807017e-02,
                     -5.1842645119807268e-02,
                     -4.7982595738891853e-02,
                     -3.1979490931262386e-02,
                     -2.9223295828573272e-02,
                     2.0540473525879417e+00,
                     -4.0028554306056982e-15,
                     -6.9951145562476089e-04,
                     2.4976240233266678e-01,
                     -2.4976240233266689e-01,
                     -1.3564222520901301e-01,
                     -1.0382299986613142e-01,
                     -6.1609103094066843e-02,
                     6.1609103094066822e-02,
                     1.0382299986613128e-01,
                     1.3564222520901303e-01,
                     3.9793805970747786e-01,
                     -3.9793805970747775e-01,
                     -3.3093711030400735e-01,
                     -2.9925907068298935e-01,
                     -3.1709752951666426e-01,
                     -2.4143986296009243e-01,
                     -2.1138807880130292e-01,
                     -1.2258597555933898e-01,
                     -9.3186919583587394e-02,
                     9.3186919583587352e-02,
                     1.2258597555933883e-01,
                     2.1138807880130264e-01,
                     2.4143986296009246e-01,
                     3.1709752951666303e-01,
                     2.9925907068298946e-01,
                     3.3093711030400719e-01,
                     -4.0829691261238581e-15,
                     5.2944315830176734e+01,
                     2.8041994412355386e-15,
                     3.0893997692497188e-01,
                     3.0893997692497210e-01,
                     3.2076449519634419e-02,
                     -4.4312272843982585e-02,
                     -2.0573071399879411e-01,
                     -2.0573071399879403e-01,
                     -4.4312272843982516e-02,
                     3.2076449519634405e-02,
                     6.1818579816101504e-01,
                     6.1818579816101504e-01,
                     1.7680506790129824e-01,
                     9.8487534610865179e-02,
                     5.1802288176394401e-02,
                     -7.0185377400300966e-02,
                     -1.5463714707615489e-01,
                     -3.6509489748646445e-01,
                     -4.7407431635068770e-01,
                     -4.7407431635068770e-01,
                     -3.6509489748646423e-01,
                     -1.5463714707615461e-01,
                     -7.0185377400300911e-02,
                     5.1802288176394526e-02,
                     9.8487534610865138e-02,
                     1.7680506790129824e-01,
                     5.5475219724409643e-02,
                     3.8741570525971813e-15,
                     6.1449219985794123e+01,
                     -6.2063822493883113e-02,
                     -6.2063822493883432e-02,
                     3.0299997507500914e-01,
                     4.0913054039110580e-01,
                     1.1156697967799378e-01,
                     1.1156697967799391e-01,
                     4.0913054039110641e-01,
                     3.0299997507501014e-01,
                     1.4543328296448935e-01,
                     1.4543328296448974e-01,
                     4.5073383572840853e-01,
                     5.4818386798303620e-01,
                     9.6107568365246854e-01,
                     1.0213763334438415e+00,
                     8.9233990560373155e-01,
                     6.7790982840294656e-01,
                     4.6608215829130242e-03,
                     4.6608215829131066e-03,
                     6.7790982840294700e-01,
                     8.9233990560373255e-01,
                     1.0213763334438437e+00,
                     9.6107568365247054e-01,
                     5.4818386798303820e-01,
                     4.5073383572841036e-01,
                     -1.0926694464024132e+01,
                     8.2596769365982100e-14,
                     -1.6341262413353828e+01,
                     4.4138588847412552e-01,
                     -4.4138588847412430e-01,
                     1.2400157825195744e-01,
                     4.5365306661180627e-01,
                     3.9762993172217476e-01,
                     -3.9762993172217503e-01,
                     -4.5365306661180715e-01,
                     -1.2400157825195748e-01,
                     9.3644829509456107e-01,
                     -9.3644829509456262e-01,
                     -4.9856751877009643e-02,
                     6.6895231508357328e-02,
                     5.8566016933660237e-01,
                     8.1050814950844985e-01,
                     1.3121994212045380e+00,
                     1.3789309401065006e+00,
                     2.4388945409002091e-01,
                     -2.4388945409002064e-01,
                     -1.3789309401065024e+00,
                     -1.3121994212045405e+00,
                     -8.1050814950845207e-01,
                     -5.8566016933660237e-01,
                     -6.6895231508356634e-02,
                     4.9856751877011288e-02,
                     1.1604284433941885e-14,
                     -4.2883667549932156e+01,
                     3.7081073854329425e-14,
                     5.5399671925113925e-01,
                     5.5399671925113836e-01,
                     -2.1376350783558574e-02,
                     1.7544652069348685e-01,
                     2.5730917965987432e-01,
                     2.5730917965987465e-01,
                     1.7544652069348698e-01,
                     -2.1376350783558571e-02,
                     3.4235719115587329e+00,
                     3.4235719115587311e+00,
                     -1.4090263701860800e-01,
                     -7.3502113187415688e-02,
                     6.7713704000015251e-02,
                     2.5302553182419463e-01,
                     5.9736487677408057e-01,
                     6.2505531159383387e-01,
                     6.7170756775031615e-01,
                     6.7170756775031648e-01,
                     6.2505531159383476e-01,
                     5.9736487677408112e-01,
                     2.5302553182419490e-01,
                     6.7713704000015265e-02,
                     -7.3502113187415577e-02,
                     -1.4090263701860770e-01,
                     -1.2778820444232183e+01,
                     2.6836027465441879e-14,
                     1.9549691287812610e+01,
                     -5.7322544947253518e-01,
                     5.7322544947253584e-01,
                     1.3446147112453297e+00,
                     9.3388295563100232e-01,
                     1.1829987547662808e-01,
                     -1.1829987547662829e-01,
                     -9.3388295563100254e-01,
                     -1.3446147112453306e+00,
                     -4.2381024604681485e-03,
                     4.2381024604675101e-03,
                     2.8239431292931236e+00,
                     2.7888700811265608e+00,
                     3.4239903101116225e+00,
                     2.6977947651927785e+00,
                     1.5547406854006427e+00,
                     5.9647783158049750e-01,
                     1.7255807744184545e-01,
                     -1.7255807744184526e-01,
                     -5.9647783158049827e-01,
                     -1.5547406854006434e+00,
                     -2.6977947651927807e+00,
                     -3.4239903101116176e+00,
                     -2.7888700811265621e+00,
                     -2.8239431292931227e+00,
                     3.2271580131232790e-14,
                     -2.5810998731834354e+02,
                     7.3411766339700258e-15,
                     -2.7738681367299871e+00,
                     -2.7738681367299880e+00,
                     -1.4730268717969874e-01,
                     1.9193318688683614e+00,
                     8.0319665180695221e-01,
                     8.0319665180695321e-01,
                     1.9193318688683618e+00,
                     -1.4730268717969813e-01,
                     2.0618387237448847e+00,
                     2.0618387237448861e+00,
                     -3.0340779319882278e+00,
                     -1.5786406928015178e+00,
                     1.4752072356445805e+00,
                     3.5599931781280953e+00,
                     4.8853919884823451e+00,
                     4.2389796663893220e+00,
                     6.6077463344223084e-02,
                     6.6077463344223222e-02,
                     4.2389796663893238e+00,
                     4.8853919884823460e+00,
                     3.5599931781280989e+00,
                     1.4752072356445818e+00,
                     -1.5786406928015175e+00,
                     -3.0340779319882274e+00,
                     -2.2952254655418667e+01,
                     4.8021605169019831e-14,
                     -2.3167356847026929e+02,
                     -4.1834803974802446e+00,
                     4.1834803974802455e+00,
                     -1.6090011489807704e-01,
                     7.8306533718787885e-01,
                     1.8123564764633626e+00,
                     -1.8123564764633626e+00,
                     -7.8306533718787885e-01,
                     1.6090011489807698e-01,
                     -8.8259411688461320e+00,
                     8.8259411688461302e+00,
                     6.7572395635899396e-01,
                     9.7500521315770625e-02,
                     -4.4944585785307034e-01,
                     4.6975646796563736e-01,
                     3.4654686215072963e+00,
                     4.9235440818664875e+00,
                     1.4525262796037648e+00,
                     -1.4525262796037643e+00,
                     -4.9235440818664884e+00,
                     -3.4654686215072967e+00,
                     -4.6975646796563830e-01,
                     4.4944585785307128e-01,
                     -9.7500521315771249e-02,
                     -6.7572395635899440e-01,
                     2.6823765829054211e-14,
                     -2.5212863247650543e+02,
                     5.1169623670909437e-16,
                     -2.1602997149717491e+00,
                     -2.1602997149717491e+00,
                     -6.7338163379221841e-02,
                     -1.1201154131599886e-01,
                     1.7449464768043175e+00,
                     1.7449464768043172e+00,
                     -1.1201154131599893e-01,
                     -6.7338163379221827e-02,
                     -1.0320753791134457e+01,
                     -1.0320753791134457e+01,
                     -5.9407363472725808e-02,
                     -5.8420218832951003e-02,
                     -4.5538576021222094e-01,
                     -3.8136052651487828e-01,
                     3.0346409006764019e-01,
                     2.1192275334091475e+00,
                     5.0872250579068616e+00,
                     5.0872250579068616e+00,
                     2.1192275334091470e+00,
                     3.0346409006763952e-01,
                     -3.8136052651487856e-01,
                     -4.5538576021222127e-01,
                     -5.8420218832951003e-02,
                     -5.9407363472726175e-02,
                     8.7202278432924771e+00,
                     -9.8632459157166315e-15,
                     -3.6291252257412879e+02)
                        .finished();

    return refPiNabla;
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_O1)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(1);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_RefPiNabla()[0];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, geometry_utilities.Tolerance1D()));

    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.6e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(4.8e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-14, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_O2)
{

    const std::string exportFolder = "Test_VEM_PCC_2D_O2";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(2);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(polygon_data.Vertices);
        vtkUtilities.Export(exportFolder + "/Polygon.vtu");

        Test_VEM_PCC_2D_Export_Dofs(polygon, reference_element_data, exportFolder);
    }

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_RefPiNabla()[1];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, std::numeric_limits<double>::epsilon()));
    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.6e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(4.8e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-14, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(3);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_RefPiNabla()[2];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, std::numeric_limits<double>::epsilon()));
    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.7e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(5.4e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-13, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-13, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}
} // namespace UnitTesting
} // namespace Polydim

#endif
