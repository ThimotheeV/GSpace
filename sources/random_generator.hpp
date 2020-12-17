#pragma once

#include <limits>
#include <random>

//for optimize draw rand
int const PRECISION = 1000000;

//Random generator,
struct rand_gen_c
{
    //delete copy operator / constructor by copy to avoid multiple instance of rand_gen in run time
    rand_gen_c(const rand_gen_c &) = delete;
    rand_gen_c(rand_gen_c &&) = delete;
    rand_gen_c &operator=(const rand_gen_c) = delete;
    rand_gen_c &operator=(rand_gen_c &&) = delete;
    ~rand_gen_c() = default;

    rand_gen_c()
    {
        Uni_int_rand = std::uniform_int_distribution<int>(0, PRECISION);

        Seed_gen.seed(0);
    }

    void put_seed(unsigned long seed)
    {
        Seed_gen.seed(seed);
    }

    int int_0_PRECISION_rand()
    {
        return Uni_int_rand(Seed_gen);
    }
    //Juste need the last bit of the value return by Seed_gen() => XXXXXY & 000001
    bool rand_bool()
    {
        return static_cast<bool>(Seed_gen() & 0b0000'0001);
    }

    int rand_nucl()
    {
        return static_cast<int>(Seed_gen() & 0b0000'0011);
    }

    std::uniform_int_distribution<int> Uni_int_rand;
    std::mt19937_64 Seed_gen;
};
