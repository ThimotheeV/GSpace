#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <array>
#include "random_generator.hpp"

TEST_CASE("random_generator_test")
{
    SECTION("construct_rand_gen_c")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        REQUIRE(rand_gen.rand_bool() == 1);
        REQUIRE(rand_gen.int_0_PRECISION_rand() == Approx(0.19576 * PRECISION).margin(0.00001 * PRECISION));
    }

    SECTION("rand_seed")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        std::mt19937_64 rand(3);

        REQUIRE(rand() == (rand_gen.Seed_gen)());

        std::uniform_int_distribution distrib(0, 1000);

        REQUIRE(distrib(rand) == distrib(rand_gen.Seed_gen));

        rand_gen.put_seed(3);
        rand.seed(3);

        REQUIRE(distrib(rand) == distrib(rand_gen.Seed_gen));
    }

    SECTION("rand_seed_pointer")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        auto rand_ptr = &rand_gen;

        std::mt19937_64 rand(3);

        REQUIRE(rand() == (rand_ptr->Seed_gen)());

        std::uniform_int_distribution distrib(0, 1000);

        REQUIRE(distrib(rand) == distrib(rand_ptr->Seed_gen));

        rand_ptr->put_seed(3);
        rand.seed(3);

        REQUIRE(distrib(rand) == distrib(rand_ptr->Seed_gen));
    }

    SECTION("rand_bool")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(7);

        std::array<int, 2> count = {0, 0};

        for (int i = 0; i < 10000; ++i)
        {
            if (rand_gen.rand_bool())
            {
                ++count.at(0);
            }
            else
            {
                ++count.at(1);
            }
        }

        REQUIRE(count.at(0) == Approx(5000).margin(50));
        REQUIRE(count.at(1) == Approx(5000).margin(50));
    }

    SECTION("rand_nucl")
    {
        rand_gen_c rand_gen;
        rand_gen.put_seed(7);

        std::array<int, 4> count = {0, 0, 0, 0};

        for (int i = 0; i < 100000; ++i)
        {
            switch (rand_gen.rand_nucl())
            {
            case 0:
                ++count.at(0);
                break;
            case 1:
                ++count.at(1);
                break;
            case 2:
                ++count.at(2);
                break;
            case 3:
                ++count.at(3);
                break;

            default:
                break;
            }
        }

        REQUIRE(count.at(0) == Approx(25000).margin(250));
        REQUIRE(count.at(1) == Approx(25000).margin(250));
        REQUIRE(count.at(2) == Approx(25000).margin(250));
        REQUIRE(count.at(3) == Approx(25000).margin(250));
    }
}