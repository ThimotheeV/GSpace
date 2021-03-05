#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "migration.hpp"

TEST_CASE("distri_test")
{
    SECTION("construct_geometric_fwd_distrib")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);

        REQUIRE(fwd_distrib[0] == fwd_distrib[1]);

        REQUIRE(fwd_distrib[0].size() == 6);
        REQUIRE(fwd_distrib[0][0] == Approx(0.80000).margin(0.00001));
        REQUIRE(fwd_distrib[0][1] / 2 == Approx(0.02974).margin(0.00001));
        REQUIRE(fwd_distrib[0][2] / 2 == Approx(0.02379).margin(0.00001));
        REQUIRE(fwd_distrib[0][3] / 2 == Approx(0.01903).margin(0.00001));
        REQUIRE(fwd_distrib[0][4] / 2 == Approx(0.01523).margin(0.00001));
        REQUIRE(fwd_distrib[0][5] / 2 == Approx(0.01218).margin(0.00001));
    }

    SECTION("construct_geometric_fwd_cumul_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c non_cumul_fwd_distrib(demo_param);
        cumul_fwd_disp_distrib_c fwd_distrib(&rand_gen, non_cumul_fwd_distrib);

        REQUIRE(fwd_distrib.Cumul_fwd_distrib[0] == fwd_distrib.Cumul_fwd_distrib[1]);

        REQUIRE(fwd_distrib.Cumul_fwd_distrib[0]->size() == 11);
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(0)) == Approx(0.80000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(1)) == Approx(0.82975 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(2)) == Approx(0.85950 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(3)) == Approx(0.88329 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(4)) == Approx(0.90709 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(5)) == Approx(0.92613 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(6)) == Approx(0.94517 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(7)) == Approx(0.96040 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(8)) == Approx(0.97563 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(9)) == Approx(0.98782 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(fwd_distrib.Cumul_fwd_distrib[0]->at(10)) == Approx(1.00000 * PRECISION).margin(0.00001 * PRECISION));
    }

    SECTION("construct_uniform_diff_x_y_fwd_distrib")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {5, 2};
        demo_param.Proba_migr = 0.5;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        fwd_disp_distrib_c fwd_distrib(demo_param);
        //Semi_distribution
        REQUIRE(fwd_distrib[0].size() == 6);
        REQUIRE(fwd_distrib[0][0] == Approx(0.5).margin(0.01));
        REQUIRE(fwd_distrib[0][1] / 2 == Approx(0.05).margin(0.001));
        REQUIRE(fwd_distrib[0][2] / 2 == Approx(0.05).margin(0.001));
        REQUIRE(fwd_distrib[0][3] / 2 == Approx(0.05).margin(0.001));
        REQUIRE(fwd_distrib[0][4] / 2 == Approx(0.05).margin(0.001));
        REQUIRE(fwd_distrib[0][5] / 2 == Approx(0.05).margin(0.001));

        REQUIRE(fwd_distrib[1].size() == 3);
        REQUIRE(fwd_distrib[1][0] == Approx(0.5).margin(0.01));
        REQUIRE(fwd_distrib[1][1] / 2 == Approx(0.125).margin(0.0001));
        REQUIRE(fwd_distrib[1][2] / 2 == Approx(0.125).margin(0.0001));
    }

    SECTION("construct_uniform_diff_x_y_fwd_cumul_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {5, 2};
        demo_param.Proba_migr = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.Lattice.size() == 121);
        REQUIRE(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->size() == 11);
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(0)) == Approx(0.80000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(1)) == Approx(0.82000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(2)) == Approx(0.84000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(3)) == Approx(0.86000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(4)) == Approx(0.88000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(5)) == Approx(0.90000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(6)) == Approx(0.92000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(7)) == Approx(0.94000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(8)) == Approx(0.96000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(9)) == Approx(0.98000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(10)) == Approx(1.00000 * PRECISION).margin(0.00001 * PRECISION));

        REQUIRE(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->size() == 5);
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(0)) == Approx(0.80000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(1)) == Approx(0.85000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(2)) == Approx(0.90000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(3)) == Approx(0.95000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(4)) == Approx(1.00000 * PRECISION).margin(0.00001 * PRECISION));
    }

    SECTION("construct_pareto_fwd_distrib")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {48, 0};
        demo_param.Proba_migr = 0.6;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::pareto(1.24609);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::pareto(1.24609);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        //Value come from IBDSim with Dispersal_Distribution = 3
        REQUIRE(fwd_distrib[0].size() == 49);
        REQUIRE(fwd_distrib[0][0] == Approx(0.4).margin(0.01));
        REQUIRE(fwd_distrib[0][1] / 2 == Approx(0.193858 / 2).margin(0.000001));
        REQUIRE(fwd_distrib[0][2] / 2 == Approx(0.0817287 / 2).margin(0.0000001));
        REQUIRE(fwd_distrib[0][3] / 2 == Approx(0.0493117 / 2).margin(0.0000001));
        REQUIRE(fwd_distrib[0][4] / 2 == Approx(0.0344561 / 2).margin(0.0000001));
        REQUIRE(fwd_distrib[0][5] / 2 == Approx(0.026092 / 2).margin(0.000001));
        //...
    }

    SECTION("construct_pareto_fwd_cumul_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {1, 10};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 49};
        demo_param.Proba_migr = 0.6;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::pareto(3.79809435);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::pareto(3.79809435);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.Lattice.size() == 10);
        REQUIRE(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->size() == 99);
        //Value come from Dispersal_Distribution = 2 in IBDSim (modifie to be a symetric distrib)
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(0)) == Approx(0.400000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(1)) == Approx(0.673308 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(2)) == Approx(0.946616 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(3)) == Approx(0.966263 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(4)) == Approx(0.985911 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(5)) == Approx(0.990123 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(6)) == Approx(0.994335 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(7)) == Approx(0.995748 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(8)) == Approx(0.997160 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(9)) == Approx(0.997765 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(10)) == Approx(0.998370 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(11)) == Approx(0.998673 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(12)) == Approx(0.998976 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(13)) == Approx(0.999145 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(14)) == Approx(0.999313 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(15)) == Approx(0.999415 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(16)) == Approx(0.999516 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(17)) == Approx(0.999581 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[1]->at(18)) == Approx(0.999646 * PRECISION).margin(0.00001 * PRECISION));
        //...
    }

    SECTION("construct_gaussian_forward_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {5, 0};
        demo_param.Proba_migr = 0.6;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::gaussian(5);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::gaussian(5);

        fwd_disp_distrib_c fw_distri(demo_param);
        //Value come gaussian.r
        REQUIRE(fw_distri[0].size() == 6);
        REQUIRE(fw_distri[0][0] == Approx(0.4).margin(0.01));
        REQUIRE(fw_distri[0][1] / 2 == Approx(0.163212).margin(0.000001));
        REQUIRE(fw_distri[0][2] / 2 == Approx(0.093264).margin(0.000001));
        REQUIRE(fw_distri[0][3] / 2 == Approx(0.034974).margin(0.000001));
        REQUIRE(fw_distri[0][4] / 2 == Approx(0.007772).margin(0.000001));
        REQUIRE(fw_distri[0][5] / 2 == Approx(0.000777).margin(0.000001));
        //...
    }

    SECTION("construct_sichel_forward_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {48, 0};
        demo_param.Lattice_size = {200, 1};
        demo_param.Proba_migr = 0.1;
        demo_param.Gamma_sichel_param = -2.15;
        demo_param.Xi_sichel_param = 20.72;
        demo_param.Omega_sichel_param = -1;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::gamma_or_sichel(demo_param, 0);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::gamma_or_sichel(demo_param, 0);

        fwd_disp_distrib_c fw_distri(demo_param);
        //Value come gaussian.r
        REQUIRE(fw_distri[0].size() == 49);
        //Value come from Dispersal_Distribution = 2 in IBDSim (modifie to be a symetric distrib)

        REQUIRE(fw_distri[0][0] == Approx(0.9).margin(0.00001));
        REQUIRE(fw_distri[0][1] == Approx(0.037113).margin(0.00001));
        REQUIRE(fw_distri[0][2] == Approx(0.025557).margin(0.00001));
        REQUIRE(fw_distri[0][3] == Approx(0.015614).margin(0.00001));
        REQUIRE(fw_distri[0][4] == Approx(0.009013).margin(0.00001));
        REQUIRE(fw_distri[0][5] == Approx(0.005116).margin(0.00001));
        REQUIRE(fw_distri[0][6] == Approx(0.002924).margin(0.00001));
        REQUIRE(fw_distri[0][7] == Approx(0.001706).margin(0.00001));
        REQUIRE(fw_distri[0][8] == Approx(0.001022).margin(0.00001));
        REQUIRE(fw_distri[0][9] == Approx(0.000630).margin(0.00001));
        //...
    }

    SECTION("draw_value")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {10, 5};
        demo_param.Proba_migr = 1;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c non_cumul_fwd_distrib(demo_param);
        cumul_fwd_disp_distrib_c fwd_distrib(&rand_gen, non_cumul_fwd_distrib);

        std::array<int, 11> result_x{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int rep = 0; rep < 100000; ++rep)
        {
            int value = fwd_distrib.draw_value(0);
            if (value < 0)
            {
                value = -value;
            }

            ++result_x[value];
        }

        //semi_distrib
        for (int dim = 1; dim <= 10; ++dim)
        {
            REQUIRE(result_x[dim] == Approx(10000).margin(500));
        }

        std::array<int, 5> result_y = std::array<int, 5>{0};
        for (int rep = 0; rep <= 100000; ++rep)
        {
            int value = fwd_distrib.draw_value(1);
            if (value < 0)
            {
                value = -value;
            }
            ++result_y[value];
        }

        for (int dim = 1; dim < 5; ++dim)
        {
            REQUIRE(result_y[dim] == Approx(20000).margin(1000));
        }
    }

    SECTION("apply_fwd_cumul_distri")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c non_cumul_fwd_distrib(demo_param);
        cumul_fwd_disp_distrib_c fwd_distrib(&rand_gen, non_cumul_fwd_distrib);

        auto result = fwd_distrib.draw_values();

        REQUIRE(result[0] == 0);
    }
}

TEST_CASE("lattice_&_extend_lattice_test")
{
    SECTION("construct_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.Lattice.size() == 121);
        REQUIRE(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->size() == 11);
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(0)) == Approx(0.80000 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(1)) == Approx(0.82975 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(2)) == Approx(0.85950 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(3)) == Approx(0.88329 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(4)) == Approx(0.90709 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(5)) == Approx(0.92613 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(6)) == Approx(0.94517 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(7)) == Approx(0.96040 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(8)) == Approx(0.97563 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(9)) == Approx(0.98782 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(lat.Cumul_fwd_distrib.Cumul_fwd_distrib[0]->at(10)) == Approx(1.00000 * PRECISION).margin(0.00001 * PRECISION));
    }

    SECTION("hash_function_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.Fwd_distrib[0].size() == 6);
        REQUIRE(lat.Fwd_distrib[1].size() == 6);

        auto key = lat.hash({1, 5});

        REQUIRE(key == 16);
    }

    SECTION("add_indivs_in_node_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        indiv_c indiv1;

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.add_indiv(&indiv1, {1, 5});

        REQUIRE(lat.node({1, 5})->Indivs_in_pop.size() == 1);
        REQUIRE(lat.node({10, 5})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({5, 1})->Indivs_in_pop.size() == 0);
        REQUIRE((indiv1.Node_lat)->Coord[0] == lat.node({1, 5})->Coord[0]);
        REQUIRE((indiv1.Node_lat)->Coord[1] == lat.node({1, 5})->Coord[1]);
    }

    SECTION("delete_indivs_in_node_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {5, 5};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        indiv_c indiv1;

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.add_indiv(&indiv1, {1, 5});

        indiv1.Node_lat->Indivs_in_pop.pop_back();
        indiv1.Node_lat = nullptr;

        REQUIRE(lat.node({1, 5})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({5, 1})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({10, 5})->Indivs_in_pop.size() == 0);
        REQUIRE(indiv1.Node_lat == nullptr);
    }

    SECTION("construct_reflecting_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 2};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE(rmap.Remap_extend_lattice.size() == 5);
        REQUIRE(rmap.Remap_extend_lattice[0].size() == 4);
    }

    SECTION("remap_absorbing_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {2, 2};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE((rmap.apply_edge_effect({-2, -2})[0]) == -1);
        REQUIRE((rmap.apply_edge_effect({-2, -2})[1]) == -1);

        REQUIRE((rmap.apply_edge_effect({-1, -1})[0]) == -1);
        REQUIRE((rmap.apply_edge_effect({-1, -1})[1]) == -1);

        REQUIRE((rmap.apply_edge_effect({2, 2})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({2, 2})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({3, 3})[0]) == -1);
        REQUIRE((rmap.apply_edge_effect({3, 3})[1]) == -1);

        REQUIRE((rmap.apply_edge_effect({4, 4})[0]) == -1);
        REQUIRE((rmap.apply_edge_effect({4, 4})[1]) == -1);
    }

    SECTION("remap_simple_reflecting_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE((rmap.apply_edge_effect({0, 0})[0]) == 0);
        REQUIRE((rmap.apply_edge_effect({0, 0})[1]) == 0);

        REQUIRE((rmap.apply_edge_effect({-1, -1})[0]) == 1);
        REQUIRE((rmap.apply_edge_effect({-1, -1})[1]) == 1);

        REQUIRE((rmap.apply_edge_effect({2, 2})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({2, 2})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({3, 3})[0]) == 1);
        REQUIRE((rmap.apply_edge_effect({3, 3})[1]) == 1);
    }

    SECTION("remap_multiple_reflecting_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {4, 4};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE((rmap.apply_edge_effect({-4, -4})[0]) == 0);
        REQUIRE((rmap.apply_edge_effect({-4, -4})[1]) == 0);

        REQUIRE((rmap.apply_edge_effect({1, 1})[0]) == 1);
        REQUIRE((rmap.apply_edge_effect({1, 1})[1]) == 1);

        REQUIRE((rmap.apply_edge_effect({5, 5})[0]) == 1);
        REQUIRE((rmap.apply_edge_effect({5, 5})[1]) == 1);
    }

    SECTION("remap_simple_circular_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE((rmap.apply_edge_effect({-1, -1})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({-1, -1})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({2, 2})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({2, 2})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({3, 3})[0]) == 0);
        REQUIRE((rmap.apply_edge_effect({3, 3})[1]) == 0);
    }

    SECTION("remap_multiple_circular_effect_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {4, 4};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        REQUIRE((rmap.apply_edge_effect({-4, -4})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({-4, -4})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({1, 1})[0]) == 1);
        REQUIRE((rmap.apply_edge_effect({1, 1})[1]) == 1);

        REQUIRE((rmap.apply_edge_effect({5, 5})[0]) == 2);
        REQUIRE((rmap.apply_edge_effect({5, 5})[1]) == 2);

        REQUIRE((rmap.apply_edge_effect({6, 6})[0]) == 0);
        REQUIRE((rmap.apply_edge_effect({6, 6})[1]) == 0);
    }

    SECTION("absorbent_apply_movment_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {20, 20};
        demo_param.Proba_migr = 0.8;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        indiv_c indiv1;

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);
        lat.add_indiv(&indiv1, {1, 5});

        //TODO changer les tests aléatoires avec attendu exact en attendu aleatoire
        auto node = apply_movment(lat, rmap, &indiv1, {-3, 6});

        REQUIRE(node->Coord[0] == 4);
        REQUIRE(node->Coord[1] == 5);

        node = apply_movment(lat, rmap, &indiv1, {12, -20});

        REQUIRE(node->Coord[0] == 4);
        REQUIRE(node->Coord[1] == 6);
    }

    SECTION("reflecting_apply_movment_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {10, 10};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        indiv_c indiv1;

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);
        lat.add_indiv(&indiv1, {1, 5});

        auto node = apply_movment(lat, rmap, &indiv1, {-3, 6});

        REQUIRE(node->Coord[0] == 2);
        REQUIRE(node->Coord[1] == 9);

        node = apply_movment(lat, rmap, &indiv1, {10, -8});

        REQUIRE(node->Coord[0] == 9);
        REQUIRE(node->Coord[1] == 3);
    }

    SECTION("circular_apply_movment_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {23, 23};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        indiv_c indiv1;

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);
        lat.add_indiv(&indiv1, {1, 5});

        auto node = apply_movment(lat, rmap, &indiv1, {-3, 6});

        REQUIRE(node->Coord[0] == 9);
        REQUIRE(node->Coord[1] == 0);

        node = apply_movment(lat, rmap, &indiv1, {23, -20});

        REQUIRE(node->Coord[0] == 2);
        REQUIRE(node->Coord[1] == 7);
    }

    SECTION("compute_migr_nb_reaching_focal_node_absorbing_effect")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::absorbing, {1, 1});

        //Without absorbing effect
        REQUIRE(result.size() == 9);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 1);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 1);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[3]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[4]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[5]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[6]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[7]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[8]) == Approx(0.01).margin(0.01));

        //With absorbing effect
        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::absorbing, {0, 0});

        REQUIRE(result.size() == 4);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[3]) == Approx(0.01).margin(0.01));

        //With absorbing effect
        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::absorbing, {2, 2});

        REQUIRE(result.size() == 4);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[3]) == Approx(0.01).margin(0.01));
    }

    SECTION("compute_migr_nb_reaching_focal_node_reflecting_effect")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::reflecting, {1, 1});

        //Without reflecting effect
        REQUIRE(result.size() == 9);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 1);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 1);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[3]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[4]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[5]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[6]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[7]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[8]) == Approx(0.01).margin(0.01));

        //With reflecting effect
        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::reflecting, {0, 0});

        REQUIRE(result.size() == 4);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.16).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.16).margin(0.01));

        REQUIRE(std::get<1>(result[3]) == Approx(0.04).margin(0.01));

        //With reflecting effect

        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::reflecting, {2, 2});

        REQUIRE(result.size() == 4);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.16).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.16).margin(0.01));

        REQUIRE(std::get<1>(result[3]) == Approx(0.04).margin(0.01));
    }

    SECTION("compute_migr_nb_reaching_focal_node_circular_effect")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::circular, {1, 1});

        //Without reflecting effect
        REQUIRE(result.size() == 9);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 1);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 1);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[3]) == Approx(0.08).margin(0.01));
        REQUIRE(std::get<1>(result[4]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[5]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[6]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[7]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[8]) == Approx(0.01).margin(0.01));

        //With reflecting effect

        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::circular, {0, 0});

        REQUIRE(result.size() == 9);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<0>(result[1])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[1])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[2])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[2])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[3])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[3])->Coord[1] == 1);
        REQUIRE(std::get<1>(result[3]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[4])->Coord[0] == 1);
        REQUIRE(std::get<0>(result[4])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[4]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[5]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[6]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[7]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[8]) == Approx(0.01).margin(0.01));

        //With reflecting effect
        result = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::circular, {2, 2});

        REQUIRE(result.size() == 9);

        REQUIRE(result.size() == 9);

        REQUIRE(std::get<0>(result[0])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[0])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[0]) == Approx(0.64).margin(0.01));

        REQUIRE(std::get<0>(result[1])->Coord[0] == 1);
        REQUIRE(std::get<0>(result[1])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[1]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[2])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[2])->Coord[1] == 1);
        REQUIRE(std::get<1>(result[2]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[3])->Coord[0] == 2);
        REQUIRE(std::get<0>(result[3])->Coord[1] == 0);
        REQUIRE(std::get<1>(result[3]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<0>(result[4])->Coord[0] == 0);
        REQUIRE(std::get<0>(result[4])->Coord[1] == 2);
        REQUIRE(std::get<1>(result[4]) == Approx(0.08).margin(0.01));

        REQUIRE(std::get<1>(result[5]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[6]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[7]) == Approx(0.01).margin(0.01));
        REQUIRE(std::get<1>(result[8]) == Approx(0.01).margin(0.01));
    }

    SECTION("cumul_normalized_bcwd_vector_extend_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto temp2 = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::reflecting, {1, 1});

        auto result = rmap.cumul_normalized_bcwd(temp2);

        REQUIRE(std::get<1>(result[0]) == Approx(0.64 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[1]) == Approx(0.72 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[2]) == Approx(0.80 * PRECISION).margin(0.01 * PRECISION));

        REQUIRE(std::get<1>(result[3]) == Approx(0.88 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[4]) == Approx(0.96 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[5]) == Approx(0.97 * PRECISION).margin(0.01 * PRECISION));

        REQUIRE(std::get<1>(result[6]) == Approx(0.98 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[7]) == Approx(0.99 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[8]) == Approx(1.00 * PRECISION).margin(0.01 * PRECISION));

        temp2 = rmap.compute_migr_nb_reaching_focal_node(edge_effect_enum::reflecting, {0, 0});

        result = rmap.cumul_normalized_bcwd(temp2);

        REQUIRE(result.size() == 4);

        REQUIRE(std::get<1>(result[0]) == Approx(0.64 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[1]) == Approx(0.80 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(result[2]) == Approx(0.96 * PRECISION).margin(0.01 * PRECISION));

        REQUIRE(std::get<1>(result[3]) == Approx(1.00000 * PRECISION).margin(0.01 * PRECISION));
    }

    // Lattice

    //One d lattice

    //  *
    //  |
    //  *

    SECTION("island_modele_2_pop_lattice")
    {
        int size = 1000000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {1, 2};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {0, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int y = 0; y <= 1; ++y)
        {
            total += lat.node({0, y})->Indivs_in_pop.size();
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(800000).margin(8000));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(200000).margin(2000));

        REQUIRE(lat.Nodes_with_lineage.size() == 2);
    }

    SECTION("island_modele_10_pop_lattice")
    {
        int size = 1000000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {1, 10};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 9};
        demo_param.Proba_migr = 0.45;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {0, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int y = 0; y <= demo_param.Disp_dist_max.at(1); ++y)
        {
            total += lat.node({0, y})->Indivs_in_pop.size();
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(550000).margin(5500));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 3})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 4})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 5})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 6})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 7})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 8})->Indivs_in_pop.size() == Approx(50000).margin(500));
        REQUIRE(lat.node({0, 9})->Indivs_in_pop.size() == Approx(50000).margin(500));

        REQUIRE(lat.Nodes_with_lineage.size() == 10);
    }

    SECTION("island_modele_11_pop_lattice")
    {
        int size = 1000000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {1, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 10};
        demo_param.Proba_migr = 0.8;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {0, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int y = 0; y <= demo_param.Disp_dist_max.at(1); ++y)
        {
            total += lat.node({0, y})->Indivs_in_pop.size();
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(200000).margin(2000));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 3})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 4})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 5})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 6})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 7})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 8})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 9})->Indivs_in_pop.size() == Approx(80000).margin(800));
        REQUIRE(lat.node({0, 10})->Indivs_in_pop.size() == Approx(80000).margin(800));

        REQUIRE(lat.Nodes_with_lineage.size() == 11);
    }
    //  *---*---*
    //  |   |   |
    //  *---*---*
    //  |   |   |
    //  *---*---*

    SECTION("migration_homogeneous_lattice")
    {
        int size = 1000000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {1, 1});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int x = 0; x <= 2; ++x)
        {
            for (int y = 0; y <= 2; ++y)
            {
                total += lat.node({x, y})->Indivs_in_pop.size();
            }
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(640000).margin(32000));
        REQUIRE(lat.node({1, 2})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({2, 0})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({2, 1})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({2, 2})->Indivs_in_pop.size() == Approx(10000).margin(500));

        REQUIRE(lat.Nodes_with_lineage.size() == 9);
    }

    SECTION("calc_back_distrib_absorbing_case_migration_heterogeneous_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto node = lat.node({0, 0});

        node->Bcwd_distrib.calc_bcwd_distrib(rmap);

        REQUIRE(node->Bcwd_distrib.Cumul_bcwd_distrib.size() == 4);

        auto bcwd_distrib_iter = node->Bcwd_distrib.Cumul_bcwd_distrib.begin();

        REQUIRE(std::get<1>(*bcwd_distrib_iter) == Approx(0.79012 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.88889 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.98765 * PRECISION).margin(0.00001 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(1.00000 * PRECISION).margin(0.00001 * PRECISION));
    }

    SECTION("simple_absorbing_migration_heterogeneous_lattice")
    {
        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {0, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int x = 0; x <= 2; ++x)
        {
            for (int y = 0; y <= 2; ++y)
            {
                total += lat.node({x, y})->Indivs_in_pop.size();
            }
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(79012).margin(1600));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(9876).margin(500));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(9876).margin(500));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(1234).margin(60));
        REQUIRE(lat.node({1, 2})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({2, 0})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({2, 1})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({2, 2})->Indivs_in_pop.size() == 0);

        REQUIRE(lat.Nodes_with_lineage.size() == 4);
    }

    SECTION("calc_back_distrib_reflecting_case_migration_heterogeneous_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.Homogeneous = false;
        extend_lattice_c rmap(lat);

        auto node = lat.node({10, 10});

        node->Bcwd_distrib.calc_bcwd_distrib(rmap);

        REQUIRE(node->Bcwd_distrib.Cumul_bcwd_distrib.size() == 4);

        auto bcwd_distrib_iter = node->Bcwd_distrib.Cumul_bcwd_distrib.begin();

        REQUIRE(std::get<1>(*bcwd_distrib_iter) == Approx(0.64 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.80 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.96 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(1.00 * PRECISION).margin(0.01 * PRECISION));
    }

    SECTION("simple_reflecting_migration_heterogeneous_lattice")
    {
        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::reflecting;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.Homogeneous = false;
        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {2, 2});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int x = 0; x <= 2; ++x)
        {
            for (int y = 0; y <= 2; ++y)
            {
                total += lat.node({x, y})->Indivs_in_pop.size();
            }
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({1, 2})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({2, 0})->Indivs_in_pop.size() == 0);
        REQUIRE(lat.node({2, 1})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({2, 2})->Indivs_in_pop.size() == Approx(64000).margin(3200));

        REQUIRE(lat.Nodes_with_lineage.size() == 4);
    }

    SECTION("calc_back_distrib_circular_case_migration_heterogeneous_lattice")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {11, 11};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.Homogeneous = false;
        extend_lattice_c rmap(lat);

        auto node = lat.node({0, 0});

        node->Bcwd_distrib.calc_bcwd_distrib(rmap);

        REQUIRE(node->Bcwd_distrib.Cumul_bcwd_distrib.size() == 9);

        auto bcwd_distrib_iter = node->Bcwd_distrib.Cumul_bcwd_distrib.begin();

        REQUIRE(std::get<1>(*bcwd_distrib_iter) == Approx(0.64 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.72 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.80 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.88 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.96 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.97 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.98 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(0.99 * PRECISION).margin(0.01 * PRECISION));
        REQUIRE(std::get<1>(*++bcwd_distrib_iter) == Approx(1.00 * PRECISION).margin(0.01 * PRECISION));
    }

    SECTION("simple_circular_migration_heterogeneous_lattice")
    {
        int size = 1000000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.Homogeneous = false;
        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {0, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int x = 0; x <= 2; ++x)
        {
            for (int y = 0; y <= 2; ++y)
            {
                total += lat.node({x, y})->Indivs_in_pop.size();
            }
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(640000).margin(32000));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({1, 2})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({2, 0})->Indivs_in_pop.size() == Approx(80000).margin(4000));
        REQUIRE(lat.node({2, 1})->Indivs_in_pop.size() == Approx(10000).margin(500));
        REQUIRE(lat.node({2, 2})->Indivs_in_pop.size() == Approx(10000).margin(500));

        REQUIRE(lat.Nodes_with_lineage.size() == 9);
    }

    SECTION("complexe_circular_migration_heterogeneous_lattice")
    {
        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.4;
        double g_geo_param = 0.3;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        lat.Homogeneous = false;
        extend_lattice_c rmap(lat);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {1, 1});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, &rmap, indiv_map);

        auto total = 0;

        for (int x = 0; x <= 2; ++x)
        {
            for (int y = 0; y <= 2; ++y)
            {
                total += lat.node({x, y})->Indivs_in_pop.size();
            }
        }
        REQUIRE(total == size);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(12000).margin(600));
        REQUIRE(lat.node({0, 2})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(12000).margin(600));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(36000).margin(1800));
        REQUIRE(lat.node({1, 2})->Indivs_in_pop.size() == Approx(12000).margin(600));
        REQUIRE(lat.node({2, 0})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({2, 1})->Indivs_in_pop.size() == Approx(12000).margin(600));
        REQUIRE(lat.node({2, 2})->Indivs_in_pop.size() == Approx(4000).margin(200));

        REQUIRE(lat.Nodes_with_lineage.size() == 9);
    }

    /**************Reading******************/

    SECTION("reading_homogeneous_deme_size_lattice")
    {
        std::vector<std::vector<int>> sbpop_matx{{5, 5, 5},
                                                 {5, 5, 5},
                                                 {5, 5, 5}};

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(sbpop_matx, rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.node({0, 0})->Subpop_size == 5);
        REQUIRE(lat.node({0, 1})->Subpop_size == 5);
        REQUIRE(lat.node({0, 2})->Subpop_size == 5);
        REQUIRE(lat.node({1, 0})->Subpop_size == 5);
        REQUIRE(lat.node({1, 1})->Subpop_size == 5);
        REQUIRE(lat.node({1, 2})->Subpop_size == 5);
        REQUIRE(lat.node({2, 0})->Subpop_size == 5);
        REQUIRE(lat.node({2, 1})->Subpop_size == 5);
        REQUIRE(lat.node({2, 2})->Subpop_size == 5);
    }

    SECTION("reading_heterogeneous_deme_size_lattice")
    {
        // Deme size
        //  3---6---9
        //  |   |   |
        //  2---5---8
        //  |   |   |
        //  1---4---7
        std::vector<std::vector<int>> sbpop_matx{{1, 2, 3},
                                                 {4, 5, 6},
                                                 {7, 8, 9}};

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {3, 3};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(sbpop_matx, rand_gen, fwd_distrib, demo_param);

        REQUIRE(lat.node({0, 0})->Subpop_size == 1);
        REQUIRE(lat.node({0, 1})->Subpop_size == 2);
        REQUIRE(lat.node({0, 2})->Subpop_size == 3);
        REQUIRE(lat.node({1, 0})->Subpop_size == 4);
        REQUIRE(lat.node({1, 1})->Subpop_size == 5);
        REQUIRE(lat.node({1, 2})->Subpop_size == 6);
        REQUIRE(lat.node({2, 0})->Subpop_size == 7);
        REQUIRE(lat.node({2, 1})->Subpop_size == 8);
        REQUIRE(lat.node({2, 2})->Subpop_size == 9);
    }

    SECTION("reading_mig_proba_matrix_lattice")
    {
        // Deme number
        //
        //  2---4
        //  |   |
        //  1---3
        //Here x,y
        std::vector<std::vector<double>> proba_mig{{0.64, 0.16, 0.16, 0.04},
                                                   {0.16, 0.64, 0.04, 0.16},
                                                   {0.16, 0.04, 0.64, 0.16},
                                                   {0.04, 0.16, 0.16, 0.64}};

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {2, 2};
        demo_param.Pop_size_per_node = 1;

        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        lattice_c lat(rand_gen, proba_mig, demo_param);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {1, 1});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, nullptr, indiv_map);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(64000).margin(3200));
    }

    SECTION("reading_mig_proba_matrix_lattice")
    {
        // Deme number
        //
        //  2---4
        //  |   |
        //  1---3
        //Here x,y
        std::vector<std::vector<double>> proba_mig{{0.64, 0.16, 0.16, 0.04},
                                                   {0.16, 0.64, 0.04, 0.16},
                                                   {0.16, 0.04, 0.64, 0.16},
                                                   {0.04, 0.16, 0.16, 0.64}};

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {2, 2};
        demo_param.Pop_size_per_node = 1;

        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        lattice_c lat(rand_gen, proba_mig, demo_param);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {1, 1});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, nullptr, indiv_map);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(64000).margin(3200));
    }

    SECTION("reading_mig_proba_matrix_lattice_bis")
    {
        // Deme number
        //
        //  2---4
        //  |   |
        //  1---3
        //Here x,y
        std::vector<std::vector<double>> proba_mig{{0.64, 0.16, 0.16, 0.04},
                                                   {0.16, 0.64, 0.04, 0.16},
                                                   {0.16, 0.04, 0.64, 0.16},
                                                   {0.04, 0.16, 0.16, 0.64}};

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::circular;
        demo_param.Lattice_size = {2, 2};
        demo_param.Pop_size_per_node = 1;

        int size = 100000;
        indiv_stock_c indiv_map(1);

        for (int i = 0; i < size; ++i)
        {
            indiv_map.add_indiv(new indiv_c);
        }

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        lattice_c lat(rand_gen, proba_mig, demo_param);

        for (auto indiv_c : indiv_map)
        {
            lat.add_indiv(indiv_c, {1, 0});
        }

        REQUIRE(lat.Nodes_with_lineage.size() == 1);

        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        lat.Nodes_with_lineage.clear();

        migration(lat, nullptr, indiv_map);

        REQUIRE(lat.node({0, 0})->Indivs_in_pop.size() == Approx(16000).margin(800));
        REQUIRE(lat.node({0, 1})->Indivs_in_pop.size() == Approx(4000).margin(200));
        REQUIRE(lat.node({1, 0})->Indivs_in_pop.size() == Approx(64000).margin(3200));
        REQUIRE(lat.node({1, 1})->Indivs_in_pop.size() == Approx(16000).margin(800));
    }
}
