#define CATCH_CONFIG_MAIN ObjectMemoryManagerTest
#include "catch.hpp"

#include "object_memory_manager.hpp"
#include "segment.hpp"

namespace unit_test
{
    struct object_memory_manager_test
    {
        std::vector<std::size_t> get_available_objs_indexes(obj_mem_manag_c<seg> &obj)
        {
            return obj.available_objs_indexes;
        }

        std::vector<seg *> get_in_use_objs(obj_mem_manag_c<seg> &obj)
        {
            return obj.in_use_objs;
        }

        std::map<seg *, std::size_t> get_hash_ptr_to_index(obj_mem_manag_c<seg> &obj)
        {
            return obj.hash_ptr_to_index;
        }
    };
} // namespace unit_test

TEST_CASE_METHOD(unit_test::object_memory_manager_test, "object_memory_manager_test")
{
    SECTION("get_new_obj")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();

        REQUIRE(pool_seg.get_index_from_obj(a) == 0);
    }

    SECTION("resize_pool")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();
        auto b = pool_seg.get_new_obj();

        REQUIRE(pool_seg.get_index_from_obj(a) == 0);
        REQUIRE(pool_seg.get_index_from_obj(b) == 1);
    }

    SECTION("double_resize_pool_of_2")
    {
        obj_mem_manag_c<seg> pool_seg;
        std::vector<seg *> seg_vec(2000);

        for (int i = 0; i < 2000; ++i)
        {
            seg_vec[i] = pool_seg.get_new_obj();
        }

        REQUIRE(pool_seg.get_index_from_obj(seg_vec[1999]) == 1999);
    }

    SECTION("free_obj_memory_manager_test")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();
        *a = seg(0, 1, 0, nullptr, nullptr);

        pool_seg.free_obj(a);

        REQUIRE(a == nullptr);

        auto available_objs_indexes = get_available_objs_indexes(pool_seg);
        auto last_index = available_objs_indexes.end();

        auto in_use_objs = get_in_use_objs(pool_seg);
        auto hash_ptr_to_index = get_hash_ptr_to_index(pool_seg);

        REQUIRE(*(--last_index) == 0);
        REQUIRE(in_use_objs[0]->R_right_brkpt == 0);
        REQUIRE(hash_ptr_to_index.size() == 1);
    }

    //WARNING : Problème de multi owner potentiel !!!!
    SECTION("get_released_index")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();

        auto b = a;

        pool_seg.free_obj(a);

        REQUIRE(a == nullptr);
        REQUIRE(pool_seg.get_index_from_obj(b) == 0);
    }

    SECTION("realocate")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();
        *a = seg(0, 1, 0, nullptr, nullptr);

        pool_seg.free_obj(a);

        auto b = pool_seg.get_new_obj();
        *b = seg(0, 2, 0, nullptr, nullptr);

        REQUIRE(pool_seg.get_index_from_obj(b) == 0);
        REQUIRE(b->R_right_brkpt == 2);
    }

    SECTION("Crash")
    {
        obj_mem_manag_c<seg> pool_seg;

        auto a = pool_seg.get_new_obj();
        auto b = pool_seg.get_new_obj();
        auto c = pool_seg.get_new_obj();

        REQUIRE(pool_seg.get_index_from_obj(c) == 2);

        auto ptr_a = a;
        pool_seg.free_obj(a);
        pool_seg.free_obj(b);

        auto d = pool_seg.get_new_obj();
        auto e = pool_seg.get_new_obj();
        auto f = pool_seg.get_new_obj();

        REQUIRE(pool_seg.get_index_from_obj(d) == 1);
        REQUIRE(pool_seg.get_index_from_obj(e) == 0);
        REQUIRE(e == ptr_a);

        REQUIRE(pool_seg.get_index_from_obj(f) == 3);
    }
}