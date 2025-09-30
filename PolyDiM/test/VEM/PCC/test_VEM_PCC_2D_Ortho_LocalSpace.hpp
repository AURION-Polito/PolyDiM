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

#ifndef __TEST_VEM_PCC_2D_Ortho_LocalSpace_H
#define __TEST_VEM_PCC_2D_Ortho_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "VEM_PCC_2D_Ortho_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_VEM_PCC_2D_Ortho_Polygon_Geometry final
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

Test_VEM_PCC_2D_Ortho_Polygon_Geometry Test_VEM_PCC_2D_Ortho_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_PCC_2D_Ortho_Polygon_Geometry result;

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

std::vector<Eigen::MatrixXd> Test_VEM_PCC_2D_Ortho_RefPiNabla()
{
    std::vector<Eigen::MatrixXd> refPiNabla(3);
    refPiNabla[0] = (Eigen::MatrixXd(3, 8) << 5.0733507212912489e-03,
                     5.0733507212912489e-03,
                     2.5606270399942312e-03,
                     3.2394650744649671e-03,
                     3.2590551741556598e-03,
                     3.2590551741556594e-03,
                     3.2394650744649667e-03,
                     2.5606270399942312e-03,
                     -1.2151417459525320e-03,
                     1.2151417459525320e-03,
                     2.9011509184616686e-03,
                     3.2201256267742084e-03,
                     1.5341164542650716e-03,
                     -1.5341164542650716e-03,
                     -3.2201256267742084e-03,
                     -2.9011509184616686e-03,
                     -4.7411388406025257e-03,
                     -4.7411388406025257e-03,
                     5.4224168558210720e-04,
                     1.5711105248917467e-03,
                     2.6277866301286732e-03,
                     2.6277866301286732e-03,
                     1.5711105248917463e-03,
                     5.4224168558210709e-04)
                        .finished();
    refPiNabla[1] = (Eigen::MatrixXd(6, 17) << -6.5524954191946664e-19,
                     -8.7712584850040846e-19,
                     1.4542663504279077e-19,
                     2.6347614929013332e-21,
                     -4.8201047287169698e-19,
                     -4.3573413878037958e-19,
                     1.5836896480014253e-19,
                     3.4341478015207322e-19,
                     -3.4384975722598988e-18,
                     1.7415811245163716e-19,
                     2.3857089606694545e-19,
                     -4.2152947346150743e-19,
                     -1.1714314243350264e-18,
                     -1.6751154445226056e-19,
                     6.3516382513233344e-19,
                     5.5052745968968608e-19,
                     7.9891000000000011e-04,
                     -1.0758457400018269e-03,
                     1.0758457400018261e-03,
                     9.8776147169394191e-04,
                     7.9128353435877407e-04,
                     2.4394510208351856e-04,
                     -2.4394510208351864e-04,
                     -7.9128353435877462e-04,
                     -9.8776147169394235e-04,
                     -8.0906419106632118e-19,
                     1.8677098546373424e-03,
                     1.9928431222609815e-03,
                     1.3041638292445546e-03,
                     -3.8621869437572026e-20,
                     -1.3041638292445553e-03,
                     -1.9928431222609819e-03,
                     -1.8677098546373433e-03,
                     1.0704993470362898e-19,
                     -1.3948536284480631e-03,
                     -1.3948536284480624e-03,
                     -1.4896540454293871e-04,
                     3.0258807851864467e-04,
                     1.0557150167337739e-03,
                     1.0557150167337737e-03,
                     3.0258807851864478e-04,
                     -1.4896540454293879e-04,
                     -5.0842392616358995e-03,
                     -4.8687972897889506e-04,
                     2.0843048534418316e-05,
                     1.3201983202650925e-03,
                     2.6338522310403343e-03,
                     1.3201983202650923e-03,
                     2.0843048534418424e-05,
                     -4.8687972897889560e-04,
                     3.1461643450756702e-05,
                     2.7434077612013803e-04,
                     2.7434077612013748e-04,
                     9.0900220583279754e-04,
                     7.7080970613562926e-04,
                     7.6842227379111952e-05,
                     7.6842227379111830e-05,
                     7.7080970613562893e-04,
                     9.0900220583279765e-04,
                     -4.4864604607001041e-04,
                     1.5313147134108753e-03,
                     1.8747261686842224e-03,
                     9.7701404408681710e-04,
                     -1.9348414442310646e-04,
                     9.7701404408681624e-04,
                     1.8747261686842224e-03,
                     1.5313147134108757e-03,
                     -3.4443637921171649e-04,
                     1.7137170101219636e-03,
                     -1.7137170101219625e-03,
                     -5.2911682335420339e-05,
                     7.2067141973420408e-04,
                     6.8320708726483318e-04,
                     -6.8320708726483297e-04,
                     -7.2067141973420376e-04,
                     5.2911682335420359e-05,
                     1.3318815329426400e-18,
                     -6.3235191144501751e-04,
                     6.5189121032484699e-04,
                     1.8938934654656907e-03,
                     2.2753928907233157e-19,
                     -1.8938934654656894e-03,
                     -6.5189121032484709e-04,
                     6.3235191144501784e-04,
                     -7.2977341582721280e-20,
                     1.2114730038724931e-03,
                     1.2114730038724928e-03,
                     -5.0313075236283184e-05,
                     2.0897286220027049e-04,
                     8.8941228265623405e-04,
                     8.8941228265623372e-04,
                     2.0897286220027047e-04,
                     -5.0313075236283157e-05,
                     4.8885244965646679e-03,
                     -2.9056713650620465e-05,
                     4.0265172892207781e-05,
                     1.0095011785356026e-03,
                     2.1082365218518164e-03,
                     1.0095011785356029e-03,
                     4.0265172892207754e-05,
                     -2.9056713650620249e-05,
                     -3.8319619505314744e-04)
                        .finished();
    refPiNabla[2] = (Eigen::MatrixXd(10, 27) << 2.6622223916816171e-20,
                     -1.8572489999184797e-19,
                     -1.4607468648073176e-19,
                     -9.1921572530143820e-20,
                     -1.5304203339900106e-19,
                     -1.4057577429385023e-19,
                     -2.2036900702227084e-19,
                     -2.9923616262636182e-19,
                     5.7459162558476005e-19,
                     7.0257320538458622e-21,
                     -4.1004101467146979e-19,
                     -3.6534841779326277e-19,
                     -3.2757776121916027e-19,
                     -2.7075726573598838e-19,
                     -2.1659946288508176e-19,
                     -2.9031158602142531e-19,
                     -3.7075588473556793e-19,
                     -3.1241944834284323e-19,
                     -4.1201063500056651e-19,
                     -4.3615016669594507e-19,
                     -6.7158106068196821e-19,
                     -7.7848194991018276e-19,
                     -6.4604720494558732e-19,
                     -6.4998851963720873e-19,
                     7.9891000000000055e-04,
                     1.6302617212973587e-20,
                     1.3145940761282824e-19,
                     -2.4447980557940692e-07,
                     2.4447980557951004e-07,
                     -3.2548442628824585e-05,
                     -3.5904266419139795e-05,
                     -5.8904276388010654e-05,
                     5.8904276388010471e-05,
                     3.5904266419138968e-05,
                     3.2548442628823304e-05,
                     6.0064261441677724e-06,
                     -6.0064261441676691e-06,
                     2.1925944715996750e-05,
                     8.7975920646303148e-06,
                     -3.1893529306597312e-05,
                     7.0747539763675615e-05,
                     -3.9286658761010358e-05,
                     1.9027648224819740e-04,
                     -2.4455693494950865e-04,
                     2.4455693494950838e-04,
                     -1.9027648224819688e-04,
                     3.9286658761008704e-05,
                     -7.0747539763673189e-05,
                     3.1893529306598044e-05,
                     -8.7975920646313973e-06,
                     -2.1925944715997574e-05,
                     4.8174511457233398e-20,
                     7.9859094059839763e-04,
                     8.8894244865193692e-20,
                     1.5760926860352683e-04,
                     1.5760926860352748e-04,
                     -4.3198045829206148e-05,
                     3.8530472580560166e-05,
                     -5.3338126079494924e-05,
                     -5.3338126079493725e-05,
                     3.8530472580560248e-05,
                     -4.3198045829206209e-05,
                     -2.5880790332569864e-04,
                     -2.5880790332569626e-04,
                     -5.7886569610563274e-05,
                     -1.6436635514489029e-04,
                     -2.1910800153542712e-06,
                     -1.0578484429011830e-04,
                     2.8755596265046552e-04,
                     -5.4096639262086793e-05,
                     1.2500805867480340e-04,
                     1.2500805867480191e-04,
                     -5.4096639262087437e-05,
                     2.8755596265046509e-04,
                     -1.0578484429011800e-04,
                     -2.1910800153535004e-06,
                     -1.6436635514489062e-04,
                     -5.7886569610563531e-05,
                     7.4034956907096855e-06,
                     3.7301995717012030e-20,
                     7.5435347652342405e-04,
                     2.7795758365547455e-04,
                     2.7795758365547460e-04,
                     4.6201177344094824e-04,
                     2.8682469185476654e-04,
                     -1.5335021679699746e-06,
                     -1.5335021679701284e-06,
                     2.8682469185476649e-04,
                     4.6201177344094846e-04,
                     -3.9369787675919604e-04,
                     -3.9369787675919599e-04,
                     1.1168429112233353e-03,
                     1.0339370198622774e-03,
                     1.1640521974447253e-03,
                     9.2034546760347020e-04,
                     4.9322941990091263e-04,
                     2.6153711926862183e-04,
                     -1.2063079051094128e-04,
                     -1.2063079051094121e-04,
                     2.6153711926862194e-04,
                     4.9322941990091263e-04,
                     9.2034546760347009e-04,
                     1.1640521974447253e-03,
                     1.0339370198622780e-03,
                     1.1168429112233356e-03,
                     -3.1096447732853331e-04,
                     4.0230418530885584e-20,
                     6.9717744875121487e-05,
                     8.4504084282766912e-04,
                     -8.4504084282766955e-04,
                     -5.4822967984404063e-05,
                     3.4566559359010566e-04,
                     3.4937891303166369e-04,
                     -3.4937891303166353e-04,
                     -3.4566559359010539e-04,
                     5.4822967984404151e-05,
                     1.6005825030153378e-03,
                     -1.6005825030153385e-03,
                     -5.3821496467362610e-04,
                     -3.6251256548347515e-04,
                     1.8123597392995269e-04,
                     5.0871128794385335e-04,
                     1.1118464598052383e-03,
                     1.2587583896351215e-03,
                     1.7956260185698555e-04,
                     -1.7956260185698539e-04,
                     -1.2587583896351208e-03,
                     -1.1118464598052376e-03,
                     -5.0871128794385324e-04,
                     -1.8123597392995301e-04,
                     3.6251256548347553e-04,
                     5.3821496467362578e-04,
                     -4.5511752776697494e-21,
                     2.5246582810062159e-05,
                     -1.5880698840889030e-19,
                     7.1145158077154532e-04,
                     7.1145158077154586e-04,
                     -1.7690210432502883e-05,
                     -3.1481525214150876e-05,
                     4.5786069127493029e-04,
                     4.5786069127493013e-04,
                     -3.1481525214151032e-05,
                     -1.7690210432502927e-05,
                     2.5104628539952445e-03,
                     2.5104628539952449e-03,
                     1.8373765866008476e-04,
                     9.1581093617710754e-05,
                     -1.5695708470204130e-04,
                     -1.7063685201638214e-04,
                     1.8490297389011772e-04,
                     5.4612884966391865e-04,
                     1.5120106182688075e-03,
                     1.5120106182688073e-03,
                     5.4612884966391757e-04,
                     1.8490297389011704e-04,
                     -1.7063685201638247e-04,
                     -1.5695708470204138e-04,
                     9.1581093617710700e-05,
                     1.8373765866008468e-04,
                     -3.2908203637855291e-04,
                     -1.5071295167956144e-19,
                     1.0163117767293034e-05,
                     -3.7194288262026604e-05,
                     3.7194288262026747e-05,
                     4.3798005322481575e-04,
                     2.7272706625718394e-04,
                     -2.5694621652437692e-05,
                     2.5694621652437533e-05,
                     -2.7272706625718394e-04,
                     -4.3798005322481581e-04,
                     3.0899680848463468e-04,
                     -3.0899680848463365e-04,
                     8.8419601578459972e-04,
                     8.9325220288655248e-04,
                     1.1166898721113000e-03,
                     8.5087938220717836e-04,
                     3.7802070564583638e-04,
                     1.8646481999565400e-05,
                     4.4037729004725263e-06,
                     -4.4037729004728871e-06,
                     -1.8646481999566091e-05,
                     -3.7802070564583648e-04,
                     -8.5087938220717923e-04,
                     -1.1166898721112993e-03,
                     -8.9325220288655248e-04,
                     -8.8419601578459961e-04,
                     -2.5190123406136561e-22,
                     -4.8557697389407762e-04,
                     1.0518985502286141e-19,
                     -5.7999980915721962e-04,
                     -5.7999980915722027e-04,
                     -3.1898389341693412e-05,
                     4.3856464108168672e-04,
                     1.4199721909416877e-04,
                     1.4199721909416907e-04,
                     4.3856464108168683e-04,
                     -3.1898389341693419e-05,
                     7.0755994226651706e-04,
                     7.0755994226651695e-04,
                     -6.8780136639887146e-04,
                     -3.5722765044624201e-04,
                     3.4564252374478621e-04,
                     8.1747728667703893e-04,
                     1.1026628381544704e-03,
                     9.1374588132686542e-04,
                     -1.0290607859268270e-04,
                     -1.0290607859268247e-04,
                     9.1374588132686651e-04,
                     1.1026628381544711e-03,
                     8.1747728667703947e-04,
                     3.4564252374478616e-04,
                     -3.5722765044624201e-04,
                     -6.8780136639887135e-04,
                     -1.5307287562601023e-04,
                     1.0454914872676658e-19,
                     -2.6511963483705561e-04,
                     -9.5115416574952615e-04,
                     9.5115416574952615e-04,
                     -3.6582175608390808e-05,
                     1.7803737241578073e-04,
                     4.1205652916398070e-04,
                     -4.1205652916398075e-04,
                     -1.7803737241578059e-04,
                     3.6582175608390835e-05,
                     -2.0066618967461967e-03,
                     2.0066618967461950e-03,
                     1.5363228578164833e-04,
                     2.2167673372654856e-05,
                     -1.0218580209753078e-04,
                     1.0680361300662944e-04,
                     7.8790734088446271e-04,
                     1.1194147022989247e-03,
                     3.3024570224780294e-04,
                     -3.3024570224780299e-04,
                     -1.1194147022989247e-03,
                     -7.8790734088446238e-04,
                     -1.0680361300662970e-04,
                     1.0218580209753104e-04,
                     -2.2167673372655060e-05,
                     -1.5363228578164868e-04,
                     -1.2784469800516495e-20,
                     -3.7553216333592127e-04,
                     -8.8290739370292710e-20,
                     -5.3089699749086356e-04,
                     -5.3089699749086366e-04,
                     -1.6548457839816702e-05,
                     -2.7527009588928758e-05,
                     4.2882329655345437e-04,
                     4.2882329655345426e-04,
                     -2.7527009588928765e-05,
                     -1.6548457839816709e-05,
                     -2.5363412130188535e-03,
                     -2.5363412130188544e-03,
                     -1.4599451491818170e-05,
                     -1.4356859169228034e-05,
                     -1.1191175517047348e-04,
                     -9.3719939453369279e-05,
                     7.4576769670737245e-05,
                     5.2080344532267730e-04,
                     1.2501934292197770e-03,
                     1.2501934292197770e-03,
                     5.2080344532267730e-04,
                     7.4576769670737353e-05,
                     -9.3719939453369374e-05,
                     -1.1191175517047359e-04,
                     -1.4356859169228043e-05,
                     -1.4599451491818243e-05,
                     6.0572154561497365e-05,
                     9.0592094864745417e-21,
                     -5.3481225467575961e-04)
                        .finished();

    return refPiNabla;
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_Ortho_O1)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Ortho_Geometry(geometry_utilities);

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
    Polydim::VEM::PCC::VEM_PCC_2D_Ortho_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(1);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_Ortho_RefPiNabla()[0];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, geometry_utilities.Tolerance1D()));

    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.6e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(4.8e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-14, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(3.5e-14, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_Ortho_O2)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Ortho_Geometry(geometry_utilities);

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
    Polydim::VEM::PCC::VEM_PCC_2D_Ortho_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(2);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_Ortho_RefPiNabla()[1];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, std::numeric_limits<double>::epsilon()));
    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.6e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(4.8e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-14, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(3.9e-13, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}

TEST(Test_VEM_PCC, Test_VEM_PCC_2D_Ortho_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_PCC_2D_Ortho_Geometry(geometry_utilities);

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
    Polydim::VEM::PCC::VEM_PCC_2D_Ortho_LocalSpace vem_local_space;

    const auto reference_element_data = vem_reference_element.Create(3);
    const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

    // Test Reference PiNabla
    const auto refPiNabla = Test_VEM_PCC_2D_Ortho_RefPiNabla()[2];
    const double relErrPiNabla = (local_space.PiNabla - refPiNabla).norm() / refPiNabla.norm();
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.5e-14, relErrPiNabla, std::numeric_limits<double>::epsilon()));
    // Test VEM performances
    Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

    const auto result =
        performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.Monomials, vem_local_space, local_space);

    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(7.0e-15, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.7e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(5.4e-14, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
    ASSERT_EQ(result.ErrorPi0km1Grad.size(), 2);
    for (unsigned int d = 0; d < 2; ++d)
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-14, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
    ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(2.6e-12, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
}
} // namespace UnitTesting
} // namespace Polydim

#endif
