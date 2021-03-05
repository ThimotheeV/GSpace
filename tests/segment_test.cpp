#define CATCH_CONFIG_MAIN SegTest
#include "catch.hpp"
#include "segment.hpp"

TEST_CASE("segment_test")
{
    SECTION("default_construct")
    {
        seg seg1;

        REQUIRE(seg1.L_left_brkpt == 0);
        REQUIRE(seg1.R_right_brkpt == 0);
        REQUIRE(seg1.U_ident_node == 0);
        REQUIRE(!seg1.Previous_seg);
        REQUIRE(!seg1.Next_seg);
    }

    SECTION("construct_partial")
    {
        seg seg1(0, 10, 1);

        REQUIRE(seg1.L_left_brkpt == 0);
        REQUIRE(seg1.R_right_brkpt == 10);
        REQUIRE(seg1.U_ident_node == 1);
        REQUIRE(!seg1.Previous_seg);
        REQUIRE(!seg1.Next_seg);
    }

    SECTION("consrtuct_complete")
    {
        seg seg1, seg2, seg3;

        seg1 = seg(0, 5, 1, nullptr, &seg2);
        seg2 = seg(6, 10, 2, &seg1, &seg3);
        seg3 = seg(10, 15, 3, &seg2, nullptr);

        REQUIRE(seg2.L_left_brkpt == 6);
        REQUIRE(seg2.R_right_brkpt == 10);
        REQUIRE(seg2.U_ident_node == 2);
        REQUIRE(seg2.Previous_seg == &seg1);
        REQUIRE(seg2.Next_seg == &seg3);
    }

    SECTION("copy")
    {
        seg seg0;
        seg seg1(0, 5, 1, &seg0, nullptr);
        seg seg2 = seg1;

        REQUIRE(seg2.L_left_brkpt == 0);
        REQUIRE(seg2.R_right_brkpt == 5);
        REQUIRE(seg2.U_ident_node == 1);
        REQUIRE(seg2.Previous_seg == &seg0);
        REQUIRE(!(seg2.Next_seg));
    }

    SECTION("clean")
    {
        seg seg0;
        seg seg1(0, 5, 1, &seg0, nullptr);

        seg1.clean();

        REQUIRE(seg1.L_left_brkpt == 0);
        REQUIRE(seg1.R_right_brkpt == 0);
        REQUIRE(seg1.U_ident_node == 0);
        REQUIRE(seg1.Previous_seg == nullptr);
        REQUIRE(seg1.Next_seg == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
    }
}

TEST_CASE("indiv_segment_test")
{
    SECTION("nbr_chr")
    {
        auto indiv = indiv_c(3, nullptr);
        REQUIRE(indiv.nbr_chr() == 3);
    }

    SECTION("update_homolog_seg")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.update_homolog_seg(chr, &seg2);

        REQUIRE(indiv.get_chr(chr, 1) != nullptr);
        REQUIRE(indiv.get_chr(chr, 0) == &seg0);
        REQUIRE(indiv.get_chr(chr, 1) == &seg2);

        REQUIRE(seg0.get_indiv() == &indiv);
        REQUIRE(seg1.get_indiv() == nullptr);
        REQUIRE(seg2.get_indiv() == &indiv);
    }

    SECTION("update_inverse_order_indiv")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.update_indiv(chr, nullptr, &seg2);

        REQUIRE(indiv.get_chr(chr, 0) == &seg2);
        REQUIRE(indiv.get_chr(chr, 1) == nullptr);
        REQUIRE(seg2.get_indiv() == &indiv);

        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
    }

    SECTION("update_indiv")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);
        seg seg3(5, 6, 3);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.update_indiv(chr, &seg2, &seg3);

        REQUIRE(indiv.get_chr(chr, 1) != nullptr);
        REQUIRE(indiv.get_chr(chr, 0) == &seg2);
        REQUIRE(indiv.get_chr(chr, 1) == &seg3);
        REQUIRE(seg2.get_indiv() == &indiv);
        REQUIRE(seg3.get_indiv() == &indiv);

        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
    }

    SECTION("swap_update_indiv")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.update_indiv(chr, &seg1, &seg0);

        REQUIRE(indiv.get_chr(chr, 1) != nullptr);
        REQUIRE(indiv.get_chr(chr, 0) == &seg1);
        REQUIRE(indiv.get_chr(chr, 1) == &seg0);
        REQUIRE(seg0.get_indiv() == &indiv);
        REQUIRE(seg1.get_indiv() == &indiv);
    }

    SECTION("delete_chr")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.untie_chr(chr, 0);
        indiv.untie_chr(chr, 1);

        REQUIRE(indiv.get_chr(chr, 0) == nullptr);
        REQUIRE(indiv.get_chr(chr, 1) == nullptr);

        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
    }

    SECTION("tricky_update_indiv")
    {
        int chr = 0;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);
        seg seg3(5, 6, 3);

        auto indiv = indiv_c(1, nullptr);
        indiv.update_indiv(chr, &seg0, &seg1);

        indiv.untie_chr(chr, 0);
        indiv.untie_chr(chr, 1);

        indiv.update_indiv(chr, &seg2, &seg3);

        REQUIRE(indiv.get_chr(chr, 1) != nullptr);
        REQUIRE(indiv.get_chr(chr, 0) == &seg2);
        REQUIRE(indiv.get_chr(chr, 1) == &seg3);
        REQUIRE(seg2.get_indiv() == &indiv);
        REQUIRE(seg3.get_indiv() == &indiv);

        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
    }

    SECTION("find_first")
    {
        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);
        seg seg3(5, 6, 3);

        seg0.Next_seg = &seg1;

        seg1.Next_seg = &seg2;
        seg1.Previous_seg = &seg0;

        seg2.Next_seg = &seg3;
        seg2.Previous_seg = &seg1;

        seg3.Previous_seg = &seg2;

        auto first_seg = find_first_seg(&seg3);

        REQUIRE(first_seg == &seg0);
        REQUIRE(seg3.Previous_seg == &seg2);
    }

    SECTION("random_coa_process_btw_seg")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(7);

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);

        std::array<seg *, 4> seg_random_vec = {&seg0, nullptr, &seg1, nullptr};
        seg_random_vec = random_coa_process_btw_seg(seg_random_vec, rand_gen);

        REQUIRE(seg_random_vec.size() == 4);
        REQUIRE(seg_random_vec[0] == &seg0);
        REQUIRE(seg_random_vec[1] == nullptr);

        seg_random_vec = {&seg0, nullptr, &seg1, nullptr};
        random_coa_process_btw_seg(seg_random_vec, rand_gen);
        REQUIRE(seg_random_vec[0] == &seg0);
        REQUIRE(seg_random_vec[1] == &seg1);
    }

    SECTION("diploid_random_coa_process_btw_seg")
    {
        rand_gen_c rand_gen;

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);
        seg seg3(5, 6, 3);

        std::array<seg *, 4> seg_random_vec = {&seg0, &seg2, &seg1, &seg3};
        seg_random_vec = random_coa_process_btw_seg(seg_random_vec, rand_gen);

        REQUIRE(seg_random_vec.size() == 4);
        REQUIRE(seg_random_vec[0] == &seg0);
        REQUIRE(seg_random_vec[1] == &seg1);
        REQUIRE(seg_random_vec[2] == &seg2);
        REQUIRE(seg_random_vec[3] == &seg3);

        seg_random_vec = {&seg0, &seg2, &seg1, &seg3};
        random_coa_process_btw_seg(seg_random_vec, rand_gen);

        REQUIRE(seg_random_vec[0] == &seg0);
        REQUIRE(seg_random_vec[1] == &seg3);
        REQUIRE(seg_random_vec[2] == &seg2);
        REQUIRE(seg_random_vec[3] == &seg1);
    }
}

TEST_CASE("indiv_stock_segment_test")
{
    SECTION("reserve_indiv_c")
    {
        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(10);

        REQUIRE(indiv_stock.capacity() == 10);
    }

    SECTION("add_indiv_c")
    {
        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(10);

        indiv_stock.add_indiv(new indiv_c());
        indiv_stock.add_indiv(new indiv_c());
        indiv_stock.add_indiv(new indiv_c());

        REQUIRE(indiv_stock.size() == 3);
        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[1]->Ident == 1);
        REQUIRE(indiv_stock[2]->Ident == 2);
    }

    SECTION("add_indiv_c without reserve")
    {
        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(2);

        indiv_stock.add_indiv(new indiv_c());
        indiv_stock.add_indiv(new indiv_c());
        indiv_stock.add_indiv(new indiv_c());

        REQUIRE(indiv_stock.size() == 3);
        REQUIRE(indiv_stock.capacity() == 6);
        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[1]->Ident == 1);
        REQUIRE(indiv_stock[2]->Ident == 2);
    }

    SECTION("clean_indiv_at_chr_c")
    {
        int chr = 0;

        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(2);

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);
        seg seg2(5, 6, 2);

        auto indiv1 = new indiv_c(1, nullptr);
        auto indiv2 = new indiv_c(1, nullptr);
        auto indiv3 = new indiv_c(1, nullptr);

        indiv1->update_indiv(chr, &seg0, &seg1);
        indiv2->update_indiv(chr, &seg2, nullptr);
        indiv3->update_indiv(chr, nullptr, nullptr);

        indiv_stock.add_indiv(indiv1);
        indiv_stock.add_indiv(indiv2);
        indiv_stock.add_indiv(indiv3);

        REQUIRE(seg0.get_indiv() == indiv_stock[0]);
        REQUIRE(seg1.get_indiv() == indiv_stock[0]);
        REQUIRE(seg2.get_indiv() == indiv_stock[1]);

        indiv_stock.clean_indiv_at_chr(0, chr);
        indiv_stock.clean_indiv_at_chr(1, chr);
        indiv_stock.clean_indiv_at_chr(2, chr);

        REQUIRE(indiv_stock.size() == 3);
        REQUIRE(indiv_stock.capacity() == 6);

        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[1]->Ident == 1);
        REQUIRE(indiv_stock[2]->Ident == 2);

        REQUIRE(indiv_stock[0]->get_chr(chr, 0) == nullptr);
        REQUIRE(indiv_stock[0]->get_chr(chr, 1) == nullptr);

        REQUIRE(indiv_stock[1]->get_chr(chr, 0) == nullptr);
        REQUIRE(indiv_stock[1]->get_chr(chr, 1) == nullptr);

        REQUIRE(indiv_stock[2]->get_chr(chr, 0) == nullptr);
        REQUIRE(indiv_stock[2]->get_chr(chr, 1) == nullptr);

        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == nullptr);
        REQUIRE(seg2.get_indiv() == nullptr);
    }

    SECTION("erase_indiv_c")
    {
        int chr = 0;

        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(2);

        seg seg0(5, 6, 3);

        auto indiv1 = new indiv_c(1, nullptr);
        auto indiv2 = new indiv_c(1, nullptr);
        auto indiv3 = new indiv_c(1, nullptr);

        indiv1->update_indiv(chr, nullptr, nullptr);
        indiv2->update_indiv(chr, &seg0, nullptr);
        indiv3->update_indiv(chr, nullptr, nullptr);

        indiv_stock.add_indiv(indiv1);
        indiv_stock.add_indiv(indiv2);
        indiv_stock.add_indiv(indiv3);

        indiv_stock.erase(1);

        REQUIRE(indiv_stock.size() == 2);
        REQUIRE(indiv_stock.capacity() == 6);
        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[1]->Ident == 1);

        REQUIRE(seg0.get_indiv() == nullptr);

        REQUIRE_THROWS(indiv_stock.erase(2));
    }

    SECTION("erase_no_consistant_indiv_c")
    {
        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(2);

        indiv_stock.add_indiv(new indiv_c());

        indiv_stock.erase(0);

        REQUIRE(indiv_stock.size() == 0);
        REQUIRE(indiv_stock.capacity() == 2);
    }

    SECTION("new_indiv_c")
    {
        int chr = 0;

        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(10);

        seg seg0(5, 6, 3);

        auto indiv = indiv_stock.new_indiv(nullptr);
        indiv->update_indiv(chr, &seg0, nullptr);

        REQUIRE(indiv_stock.size() == 1);
        REQUIRE(indiv_stock.capacity() == 10);
        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[0] == indiv);
        REQUIRE(indiv_stock[0]->get_chr(chr, 0) == &seg0);
        REQUIRE(seg0.get_indiv() == indiv_stock[0]);
    }

    SECTION("recycle_indiv_c")
    {
        int chr = 0;

        indiv_stock_c indiv_stock(1);
        indiv_stock.reserve(10);

        seg seg0(5, 6, 0);
        seg seg1(5, 6, 1);

        auto indiv = indiv_stock.new_indiv(nullptr);
        indiv->update_indiv(chr, &seg0, nullptr);

        indiv_stock.erase(0);

        auto indiv1 = indiv_stock.new_indiv(nullptr);
        indiv1->update_indiv(chr, &seg1, nullptr);

        REQUIRE(indiv_stock.size() == 1);
        REQUIRE(indiv_stock.capacity() == 10);
        REQUIRE(indiv_stock[0]->Ident == 0);
        REQUIRE(indiv_stock[0]->get_chr(chr, 0) == &seg1);
        REQUIRE(seg0.get_indiv() == nullptr);
        REQUIRE(seg1.get_indiv() == indiv_stock[0]);
    }
}

TEST_CASE("indiv_stock_intergrate_segment_test")
{
    SECTION("ini")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 10;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 10;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_stock_c indiv_stock(samp_param.Chr_nbr);
        indiv_stock.ini(lat, samp_param);

        REQUIRE(indiv_stock.size() == static_cast<std::size_t>(samp_param.n_total_sample_size));
        REQUIRE(indiv_stock.capacity() == 10 * static_cast<std::size_t>(samp_param.n_total_sample_size));

        for (int indiv = 0; indiv < samp_param.n_total_sample_size; ++indiv)
        {
            REQUIRE(indiv_stock[indiv]->Ident == indiv);
            REQUIRE(indiv_stock[indiv]->Node_lat == lat.node({0, 0}));
            REQUIRE(indiv_stock[indiv]->nbr_chr() ==  static_cast<std::size_t>(samp_param.Chr_nbr));

            for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
            {
                REQUIRE(indiv_stock[indiv]->get_chr(chr, 0) == nullptr);
                REQUIRE(indiv_stock[indiv]->get_chr(chr, 1) == nullptr);
            }
        }
    }
}
