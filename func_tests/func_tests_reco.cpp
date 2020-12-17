#define CATCH_CONFIG_MAIN FunctionnalTestRecombinaison
#include "catch.hpp"

#include "simulator.hpp"
#include <iostream>

namespace unit_test
{
    struct FunctionnalTest
    {
        // FunctionnalTest() {}
        // ~FunctionnalTest() {}

        static double calc_phiA(double N, double m)
        {
            double result = -(std::pow(m, 2) - 2 * m + 1) / ((2 * N - 1) * std::pow(m, 2) - 2 * (2 * N - 1) * m - 1);
            return result;
        }

        static double calc_phiAB(double N, double m, double r)
        {
            double result = -((16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 14) - 14 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 13) + (1440 * std::pow(N, 6) - 7208 * std::pow(N, 5) + 14436 * std::pow(N, 4) - 14822 * std::pow(N, 3) + 8241 * std::pow(N, 2) - 2360 * N + 273) * std::pow(m, 12) - 4 * (1408 * std::pow(N, 6) - 7064 * std::pow(N, 5) + 14188 * std::pow(N, 4) - 14618 * std::pow(N, 3) + 8161 * std::pow(N, 2) - 2348 * N + 273) * std::pow(m, 11) + (14928 * std::pow(N, 6) - 75240 * std::pow(N, 5) + 151888 * std::pow(N, 4) - 157394 * std::pow(N, 3) + 88445 * std::pow(N, 2) - 25630 * N + 3003) * std::pow(m, 10) - 2 * (14096 * std::pow(N, 6) - 71720 * std::pow(N, 5) + 146080 * std::pow(N, 4) - 152754 * std::pow(N, 3) + 86661 * std::pow(N, 2) - 25366 * N + 3003) * std::pow(m, 9) + (38720 * std::pow(N, 6) - 200712 * std::pow(N, 5) + 415192 * std::pow(N, 4) - 440304 * std::pow(N, 3) + 253203 * std::pow(N, 2) - 75108 * N + 9009) * std::pow(m, 8) - 8 * (4832 * std::pow(N, 6) - 25944 * std::pow(N, 5) + 55096 * std::pow(N, 4) - 59700 * std::pow(N, 3) + 34989 * std::pow(N, 2) - 10560 * N + 1287) * std::pow(m, 7) + (27456 * std::pow(N, 6) - 157256 * std::pow(N, 5) + 348880 * std::pow(N, 4) - 390600 * std::pow(N, 3) + 235053 * std::pow(N, 2) - 72534 * N + 9009) * std::pow(m, 6) - 2 * (6592 * std::pow(N, 6) - 42520 * std::pow(N, 5) + 101424 * std::pow(N, 4) - 119320 * std::pow(N, 3) + 74495 * std::pow(N, 2) - 23650 * N + 3003) * std::pow(m, 5) + (3840 * std::pow(N, 6) - 31120 * std::pow(N, 5) + 83940 * std::pow(N, 4) - 106534 * std::pow(N, 3) + 70031 * std::pow(N, 2) - 23056 * N + 3003) * std::pow(m, 4) - 4 * (128 * std::pow(N, 6) - 1728 * std::pow(N, 5) + 5828 * std::pow(N, 4) - 8342 * std::pow(N, 3) + 5901 * std::pow(N, 2) - 2036 * N + 273) * std::pow(m, 3) - ((16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 14) - 14 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 13) + (1440 * std::pow(N, 6) - 7208 * std::pow(N, 5) + 14436 * std::pow(N, 4) - 14822 * std::pow(N, 3) + 8241 * std::pow(N, 2) - 2360 * N + 273) * std::pow(m, 12) - 4 * (1408 * std::pow(N, 6) - 7064 * std::pow(N, 5) + 14188 * std::pow(N, 4) - 14618 * std::pow(N, 3) + 8161 * std::pow(N, 2) - 2348 * N + 273) * std::pow(m, 11) + (14944 * std::pow(N, 6) - 75272 * std::pow(N, 5) + 151924 * std::pow(N, 4) - 157430 * std::pow(N, 3) + 88465 * std::pow(N, 2) - 25634 * N + 3003) * std::pow(m, 10) - 2 * (14176 * std::pow(N, 6) - 71880 * std::pow(N, 5) + 146260 * std::pow(N, 4) - 152934 * std::pow(N, 3) + 86761 * std::pow(N, 2) - 25386 * N + 3003) * std::pow(m, 9) + (39424 * std::pow(N, 6) - 202144 * std::pow(N, 5) + 416844 * std::pow(N, 4) - 441956 * std::pow(N, 3) + 254111 * std::pow(N, 2) - 75288 * N + 9009) * std::pow(m, 8) - 8 * (5056 * std::pow(N, 6) - 26416 * std::pow(N, 5) + 55668 * std::pow(N, 4) - 60272 * std::pow(N, 3) + 35297 * std::pow(N, 2) - 10620 * N + 1287) * std::pow(m, 7) + 7 * (4336 * std::pow(N, 6) - 23392 * std::pow(N, 5) + 51048 * std::pow(N, 4) - 57008 * std::pow(N, 3) + 34211 * std::pow(N, 2) - 10482 * N + 1287) * std::pow(m, 6) - 2 * (8112 * std::pow(N, 6) - 46304 * std::pow(N, 5) + 106856 * std::pow(N, 4) - 124752 * std::pow(N, 3) + 77239 * std::pow(N, 2) - 24154 * N + 3003) * std::pow(m, 5) - 16 * std::pow(N, 5) + (5856 * std::pow(N, 6) - 37160 * std::pow(N, 5) + 93740 * std::pow(N, 4) - 116334 * std::pow(N, 3) + 74791 * std::pow(N, 2) - 23896 * N + 3003) * std::pow(m, 4) + 68 * std::pow(N, 4) - 4 * (320 * std::pow(N, 6) - 2536 * std::pow(N, 5) + 7356 * std::pow(N, 4) - 9870 * std::pow(N, 3) + 6613 * std::pow(N, 2) - 2156 * N + 273) * std::pow(m, 3) - 104 * std::pow(N, 3) + (128 * std::pow(N, 6) - 1800 * std::pow(N, 5) + 6388 * std::pow(N, 4) - 9334 * std::pow(N, 3) + 6495 * std::pow(N, 2) - 2150 * N + 273) * std::pow(m, 2) + 73 * std::pow(N, 2) + 2 * (104 * std::pow(N, 5) - 452 * std::pow(N, 4) + 702 * std::pow(N, 3) - 499 * std::pow(N, 2) + 166 * N - 21) * m - 24 * N + 3) * std::pow(r, 3) - 36 * std::pow(N, 3) - (704 * std::pow(N, 5) - 3872 * std::pow(N, 4) + 6818 * std::pow(N, 3) - 5371 * std::pow(N, 2) + 1970 * N - 273) * std::pow(m, 2) + (3 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 14) - 42 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 13) + 3 * (1440 * std::pow(N, 6) - 7208 * std::pow(N, 5) + 14436 * std::pow(N, 4) - 14822 * std::pow(N, 3) + 8241 * std::pow(N, 2) - 2360 * N + 273) * std::pow(m, 12) - 12 * (1408 * std::pow(N, 6) - 7064 * std::pow(N, 5) + 14188 * std::pow(N, 4) - 14618 * std::pow(N, 3) + 8161 * std::pow(N, 2) - 2348 * N + 273) * std::pow(m, 11) + (44816 * std::pow(N, 6) - 225792 * std::pow(N, 5) + 455760 * std::pow(N, 4) - 472280 * std::pow(N, 3) + 265387 * std::pow(N, 2) - 76900 * N + 9009) * std::pow(m, 10) - 2 * (42448 * std::pow(N, 6) - 215520 * std::pow(N, 5) + 438720 * std::pow(N, 4) - 458752 * std::pow(N, 3) + 260243 * std::pow(N, 2) - 76148 * N + 9009) * std::pow(m, 9) + (117568 * std::pow(N, 6) - 605368 * std::pow(N, 5) + 1249980 * std::pow(N, 4) - 1325402 * std::pow(N, 3) + 761969 * std::pow(N, 2) - 225774 * N + 27027) * std::pow(m, 8) - 8 * (14944 * std::pow(N, 6) - 78904 * std::pow(N, 5) + 166812 * std::pow(N, 4) - 180650 * std::pow(N, 3) + 105767 * std::pow(N, 2) - 31830 * N + 3861) * std::pow(m, 7) + (88160 * std::pow(N, 6) - 486640 * std::pow(N, 5) + 1069148 * std::pow(N, 4) - 1194620 * std::pow(N, 3) + 716639 * std::pow(N, 2) - 219702 * N + 27027) * std::pow(m, 6) - 2 * (22816 * std::pow(N, 6) - 136336 * std::pow(N, 5) + 318708 * std::pow(N, 4) - 372548 * std::pow(N, 3) + 230597 * std::pow(N, 2) - 72210 * N + 9009) * std::pow(m, 5) - 32 * std::pow(N, 5) + (15552 * std::pow(N, 6) - 107552 * std::pow(N, 5) + 277800 * std::pow(N, 4) - 345782 * std::pow(N, 3) + 222413 * std::pow(N, 2) - 71268 * N + 9009) * std::pow(m, 4) + 176 * std::pow(N, 4) - 4 * (768 * std::pow(N, 6) - 7104 * std::pow(N, 5) + 21520 * std::pow(N, 4) - 29086 * std::pow(N, 3) + 19543 * std::pow(N, 2) - 6408 * N + 819) * std::pow(m, 3) - 286 * std::pow(N, 3) + (256 * std::pow(N, 6) - 4720 * std::pow(N, 5) + 18228 * std::pow(N, 4) - 27104 * std::pow(N, 3) + 19013 * std::pow(N, 2) - 6360 * N + 819) * std::pow(m, 2) + 207 * std::pow(N, 2) + 2 * (240 * std::pow(N, 5) - 1236 * std::pow(N, 4) + 1992 * std::pow(N, 3) - 1441 * std::pow(N, 2) + 488 * N - 63) * m - 70 * N + 9) * std::pow(r, 2) + 45 * std::pow(N, 2) - 2 * (144 * std::pow(N, 4) - 394 * std::pow(N, 3) + 367 * std::pow(N, 2) - 146 * N + 21) * m - (3 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 14) - 42 * (16 * std::pow(N, 6) - 80 * std::pow(N, 5) + 160 * std::pow(N, 4) - 164 * std::pow(N, 3) + 91 * std::pow(N, 2) - 26 * N + 3) * std::pow(m, 13) + 3 * (1440 * std::pow(N, 6) - 7208 * std::pow(N, 5) + 14436 * std::pow(N, 4) - 14822 * std::pow(N, 3) + 8241 * std::pow(N, 2) - 2360 * N + 273) * std::pow(m, 12) - 12 * (1408 * std::pow(N, 6) - 7064 * std::pow(N, 5) + 14188 * std::pow(N, 4) - 14618 * std::pow(N, 3) + 8161 * std::pow(N, 2) - 2348 * N + 273) * std::pow(m, 11) + (44800 * std::pow(N, 6) - 225760 * std::pow(N, 5) + 455724 * std::pow(N, 4) - 472244 * std::pow(N, 3) + 265367 * std::pow(N, 2) - 76896 * N + 9009) * std::pow(m, 10) - 2 * (42368 * std::pow(N, 6) - 215360 * std::pow(N, 5) + 438540 * std::pow(N, 4) - 458572 * std::pow(N, 3) + 260143 * std::pow(N, 2) - 76128 * N + 9009) * std::pow(m, 9) + (116864 * std::pow(N, 6) - 603936 * std::pow(N, 5) + 1248328 * std::pow(N, 4) - 1323750 * std::pow(N, 3) + 761061 * std::pow(N, 2) - 225594 * N + 27027) * std::pow(m, 8) - 8 * (14720 * std::pow(N, 6) - 78432 * std::pow(N, 5) + 166240 * std::pow(N, 4) - 180078 * std::pow(N, 3) + 105459 * std::pow(N, 2) - 31770 * N + 3861) * std::pow(m, 7) + (85264 * std::pow(N, 6) - 480144 * std::pow(N, 5) + 1060688 * std::pow(N, 4) - 1186164 * std::pow(N, 3) + 712215 * std::pow(N, 2) - 218862 * N + 27027) * std::pow(m, 6) - 2 * (21296 * std::pow(N, 6) - 132528 * std::pow(N, 5) + 313264 * std::pow(N, 4) - 367116 * std::pow(N, 3) + 227853 * std::pow(N, 2) - 71706 * N + 9009) * std::pow(m, 5) + (13536 * std::pow(N, 6) - 101384 * std::pow(N, 5) + 267940 * std::pow(N, 4) - 335982 * std::pow(N, 3) + 217653 * std::pow(N, 2) - 70428 * N + 9009) * std::pow(m, 4) + 104 * std::pow(N, 4) - 4 * (576 * std::pow(N, 6) - 6248 * std::pow(N, 5) + 19972 * std::pow(N, 4) - 27558 * std::pow(N, 3) + 18831 * std::pow(N, 2) - 6288 * N + 819) * std::pow(m, 3) - 218 * std::pow(N, 3) + (128 * std::pow(N, 6) - 3456 * std::pow(N, 5) + 15652 * std::pow(N, 4) - 24588 * std::pow(N, 3) + 17889 * std::pow(N, 2) - 6180 * N + 819) * std::pow(m, 2) + 179 * std::pow(N, 2) + 2 * (96 * std::pow(N, 5) - 916 * std::pow(N, 4) + 1684 * std::pow(N, 3) - 1309 * std::pow(N, 2) + 468 * N - 63) * m - 66 * N + 9) * r - 20 * N + 3) / ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29824 * std::pow(N, 7) - 165248 * std::pow(N, 6) + 378696 * std::pow(N, 5) - 466348 * std::pow(N, 4) + 334102 * std::pow(N, 3) - 139653 * std::pow(N, 2) + 31630 * N - 3003) * std::pow(m, 11) - 2 * (28032 * std::pow(N, 7) - 156736 * std::pow(N, 6) + 362280 * std::pow(N, 5) - 449948 * std::pow(N, 4) + 325166 * std::pow(N, 3) - 137133 * std::pow(N, 2) + 31342 * N - 3003) * std::pow(m, 9) + (76032 * std::pow(N, 7) - 433088 * std::pow(N, 6) + 1016944 * std::pow(N, 5) - 1281244 * std::pow(N, 4) + 938600 * std::pow(N, 3) - 401091 * std::pow(N, 2) + 92856 * N - 9009) * std::pow(m, 8) - 8 * (9216 * std::pow(N, 7) - 54464 * std::pow(N, 6) + 131584 * std::pow(N, 5) - 169780 * std::pow(N, 4) + 127028 * std::pow(N, 3) - 55341 * std::pow(N, 2) + 13044 * N - 1287) * std::pow(m, 7) + (49152 * std::pow(N, 7) - 312576 * std::pow(N, 6) + 794936 * std::pow(N, 5) - 1066988 * std::pow(N, 4) + 824738 * std::pow(N, 3) - 369537 * std::pow(N, 2) + 89292 * N - 9009) * std::pow(m, 6) - 2 * (10240 * std::pow(N, 7) - 76032 * std::pow(N, 6) + 212520 * std::pow(N, 5) - 304676 * std::pow(N, 4) + 247654 * std::pow(N, 3) - 115579 * std::pow(N, 2) + 28900 * N - 3003) * std::pow(m, 5) + (4096 * std::pow(N, 7) - 45056 * std::pow(N, 6) + 151696 * std::pow(N, 5) - 243524 * std::pow(N, 4) + 214156 * std::pow(N, 3) - 106063 * std::pow(N, 2) + 27802 * N - 3003) * std::pow(m, 4) + 4 * (1536 * std::pow(N, 6) - 8032 * std::pow(N, 5) + 15884 * std::pow(N, 4) - 15864 * std::pow(N, 3) + 8581 * std::pow(N, 2) - 2402 * N + 273) * std::pow(m, 3) - ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29888 * std::pow(N, 7) - 165472 * std::pow(N, 6) + 379032 * std::pow(N, 5) - 466628 * std::pow(N, 4) + 334238 * std::pow(N, 3) - 139689 * std::pow(N, 2) + 31634 * N - 3003) * std::pow(m, 10) - 2 * (28352 * std::pow(N, 7) - 157856 * std::pow(N, 6) + 363960 * std::pow(N, 5) - 451348 * std::pow(N, 4) + 325846 * std::pow(N, 3) - 137313 * std::pow(N, 2) + 31362 * N - 3003) * std::pow(m, 9) + (78848 * std::pow(N, 7) - 442976 * std::pow(N, 6) + 1031824 * std::pow(N, 5) - 1293684 * std::pow(N, 4) + 944664 * std::pow(N, 3) - 402703 * std::pow(N, 2) + 93036 * N - 9009) * std::pow(m, 8) - 8 * (10112 * std::pow(N, 7) - 57632 * std::pow(N, 6) + 136384 * std::pow(N, 5) - 173820 * std::pow(N, 4) + 129012 * std::pow(N, 3) - 55873 * std::pow(N, 2) + 13104 * N - 1287) * std::pow(m, 7) + 7 * (8672 * std::pow(N, 7) - 50576 * std::pow(N, 6) + 122656 * std::pow(N, 5) - 160176 * std::pow(N, 4) + 121674 * std::pow(N, 3) - 53839 * std::pow(N, 2) + 12876 * N - 1287) * std::pow(m, 6) - 2 * (16224 * std::pow(N, 7) - 98256 * std::pow(N, 6) + 247584 * std::pow(N, 5) - 335248 * std::pow(N, 4) + 263186 * std::pow(N, 3) - 119891 * std::pow(N, 2) + 29404 * N - 3003) * std::pow(m, 5) + 16 * std::pow(N, 5) + (11712 * std::pow(N, 7) - 75696 * std::pow(N, 6) + 202800 * std::pow(N, 5) - 290008 * std::pow(N, 4) + 238616 * std::pow(N, 3) - 113063 * std::pow(N, 2) + 28642 * N - 3003) * std::pow(m, 4) - 68 * std::pow(N, 4) - 4 * (640 * std::pow(N, 7) - 4688 * std::pow(N, 6) + 13936 * std::pow(N, 5) - 21688 * std::pow(N, 4) + 19100 * std::pow(N, 3) - 9549 * std::pow(N, 2) + 2522 * N - 273) * std::pow(m, 3) + 104 * std::pow(N, 3) + (256 * std::pow(N, 7) - 2560 * std::pow(N, 6) + 9272 * std::pow(N, 5) - 16580 * std::pow(N, 4) + 16162 * std::pow(N, 3) - 8703 * std::pow(N, 2) + 2426 * N - 273) * std::pow(m, 2) - 73 * std::pow(N, 2) + 2 * (64 * std::pow(N, 6) - 376 * std::pow(N, 5) + 868 * std::pow(N, 4) - 994 * std::pow(N, 3) + 595 * std::pow(N, 2) - 178 * N + 21) * m + 24 * N - 3) * std::pow(r, 3) + 36 * std::pow(N, 3) + (3008 * std::pow(N, 5) - 9456 * std::pow(N, 4) + 11790 * std::pow(N, 3) - 7307 * std::pow(N, 2) + 2246 * N - 273) * std::pow(m, 2) + (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89632 * std::pow(N, 7) - 496320 * std::pow(N, 6) + 1136968 * std::pow(N, 5) - 1399780 * std::pow(N, 4) + 1002660 * std::pow(N, 3) - 419051 * std::pow(N, 2) + 94900 * N - 9009) * std::pow(m, 10) - 2 * (84896 * std::pow(N, 7) - 473088 * std::pow(N, 6) + 1091240 * std::pow(N, 5) - 1353524 * std::pow(N, 4) + 977268 * std::pow(N, 3) - 411859 * std::pow(N, 2) + 94076 * N - 9009) * std::pow(m, 9) + (235136 * std::pow(N, 7) - 1324688 * std::pow(N, 6) + 3089800 * std::pow(N, 5) - 3876432 * std::pow(N, 4) + 2831586 * std::pow(N, 3) - 1207393 * std::pow(N, 2) + 279018 * N - 27027) * std::pow(m, 8) - 8 * (29888 * std::pow(N, 7) - 171536 * std::pow(N, 6) + 407320 * std::pow(N, 5) - 519960 * std::pow(N, 4) + 386250 * std::pow(N, 3) - 167383 * std::pow(N, 2) + 39282 * N - 3861) * std::pow(m, 7) + (176320 * std::pow(N, 7) - 1044240 * std::pow(N, 6) + 2551416 * std::pow(N, 5) - 3343560 * std::pow(N, 4) + 2544490 * std::pow(N, 3) - 1127371 * std::pow(N, 2) + 269976 * N - 27027) * std::pow(m, 6) - 2 * (45632 * std::pow(N, 7) - 285104 * std::pow(N, 6) + 729256 * std::pow(N, 5) - 994392 * std::pow(N, 4) + 783438 * std::pow(N, 3) - 357769 * std::pow(N, 2) + 87960 * N - 9009) * std::pow(m, 5) + 32 * std::pow(N, 5) + (31104 * std::pow(N, 7) - 213440 * std::pow(N, 6) + 588488 * std::pow(N, 5) - 852736 * std::pow(N, 4) + 706248 * std::pow(N, 3) - 336109 * std::pow(N, 2) + 85506 * N - 9009) * std::pow(m, 4) - 176 * std::pow(N, 4) - 4 * (1536 * std::pow(N, 7) - 12576 * std::pow(N, 6) + 39448 * std::pow(N, 5) - 62896 * std::pow(N, 4) + 56036 * std::pow(N, 3) - 28223 * std::pow(N, 2) + 7506 * N - 819) * std::pow(m, 3) + 286 * std::pow(N, 3) + (512 * std::pow(N, 7) - 6272 * std::pow(N, 6) + 25168 * std::pow(N, 5) - 47052 * std::pow(N, 4) + 46788 * std::pow(N, 3) - 25501 * std::pow(N, 2) + 7188 * N - 819) * std::pow(m, 2) - 207 * std::pow(N, 2) + 2 * (128 * std::pow(N, 6) - 944 * std::pow(N, 5) + 2380 * std::pow(N, 4) - 2820 * std::pow(N, 3) + 1721 * std::pow(N, 2) - 524 * N + 63) * m + 70 * N - 9) * std::pow(r, 2) - 45 * std::pow(N, 2) + 2 * (288 * std::pow(N, 4) - 574 * std::pow(N, 3) + 447 * std::pow(N, 2) - 158 * N + 21) * m - (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89568 * std::pow(N, 7) - 496096 * std::pow(N, 6) + 1136632 * std::pow(N, 5) - 1399500 * std::pow(N, 4) + 1002524 * std::pow(N, 3) - 419015 * std::pow(N, 2) + 94896 * N - 9009) * std::pow(m, 10) - 2 * (84576 * std::pow(N, 7) - 471968 * std::pow(N, 6) + 1089560 * std::pow(N, 5) - 1352124 * std::pow(N, 4) + 976588 * std::pow(N, 3) - 411679 * std::pow(N, 2) + 94056 * N - 9009) * std::pow(m, 9) + (232320 * std::pow(N, 7) - 1314800 * std::pow(N, 6) + 3074920 * std::pow(N, 5) - 3863992 * std::pow(N, 4) + 2825522 * std::pow(N, 3) - 1205781 * std::pow(N, 2) + 278838 * N - 27027) * std::pow(m, 8) - 8 * (28992 * std::pow(N, 7) - 168368 * std::pow(N, 6) + 402520 * std::pow(N, 5) - 515920 * std::pow(N, 4) + 384266 * std::pow(N, 3) - 166851 * std::pow(N, 2) + 39222 * N - 3861) * std::pow(m, 7) + (164736 * std::pow(N, 7) - 1002736 * std::pow(N, 6) + 2487736 * std::pow(N, 5) - 3289312 * std::pow(N, 4) + 2517510 * std::pow(N, 3) - 1120035 * std::pow(N, 2) + 269136 * N - 27027) * std::pow(m, 6) - 2 * (39552 * std::pow(N, 7) - 262736 * std::pow(N, 6) + 694120 * std::pow(N, 5) - 963808 * std::pow(N, 4) + 767906 * std::pow(N, 3) - 353457 * std::pow(N, 2) + 87456 * N - 9009) * std::pow(m, 5) + (23040 * std::pow(N, 7) - 182112 * std::pow(N, 6) + 537032 * std::pow(N, 5) - 806192 * std::pow(N, 4) + 681788 * std::pow(N, 3) - 329109 * std::pow(N, 2) + 84666 * N - 9009) * std::pow(m, 4) - 104 * std::pow(N, 4) - 12 * (256 * std::pow(N, 7) - 3072 * std::pow(N, 6) + 11144 * std::pow(N, 5) - 19024 * std::pow(N, 4) + 17600 * std::pow(N, 3) - 9085 * std::pow(N, 2) + 2462 * N - 273) * std::pow(m, 3) + 218 * std::pow(N, 3) - (3200 * std::pow(N, 6) - 18592 * std::pow(N, 5) + 39868 * std::pow(N, 4) - 42416 * std::pow(N, 3) + 24105 * std::pow(N, 2) - 7008 * N + 819) * std::pow(m, 2) - 179 * std::pow(N, 2) - 2 * (512 * std::pow(N, 5) - 1788 * std::pow(N, 4) + 2400 * std::pow(N, 3) - 1573 * std::pow(N, 2) + 504 * N - 63) * m + 66 * N - 9) * r + 20 * N - 3);
            return result;
        }

        static double calc_gammaAB(double N, double m, double r)
        {
            double result = ((16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 14) - 14 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 13) + 91 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 12) - 364 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 11) + (15992 * std::pow(N, 5) - 64020 * std::pow(N, 4) + 96058 * std::pow(N, 3) - 68051 * std::pow(N, 2) + 23020 * N - 3003) * std::pow(m, 10) - 2 * (15896 * std::pow(N, 5) - 63844 * std::pow(N, 4) + 95906 * std::pow(N, 3) - 67983 * std::pow(N, 2) + 23008 * N - 3003) * std::pow(m, 9) + (46960 * std::pow(N, 5) - 190188 * std::pow(N, 4) + 286556 * std::pow(N, 3) - 203433 * std::pow(N, 2) + 68934 * N - 9009) * std::pow(m, 8) - 8 * (6496 * std::pow(N, 5) - 26772 * std::pow(N, 4) + 40592 * std::pow(N, 3) - 28911 * std::pow(N, 2) + 9822 * N - 1287) * std::pow(m, 7) + (42792 * std::pow(N, 5) - 182284 * std::pow(N, 4) + 279692 * std::pow(N, 3) - 200466 * std::pow(N, 2) + 68439 * N - 9009) * std::pow(m, 6) - 2 * (12792 * std::pow(N, 5) - 57860 * std::pow(N, 4) + 90692 * std::pow(N, 3) - 65758 * std::pow(N, 2) + 22645 * N - 3003) * std::pow(m, 5) + (10544 * std::pow(N, 5) - 53204 * std::pow(N, 4) + 86576 * std::pow(N, 3) - 64078 * std::pow(N, 2) + 22393 * N - 3003) * std::pow(m, 4) - 4 * (672 * std::pow(N, 5) - 4188 * std::pow(N, 4) + 7288 * std::pow(N, 3) - 5594 * std::pow(N, 2) + 2003 * N - 273) * std::pow(m, 3) - ((16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 14) - 14 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 13) + 91 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 12) - 364 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 11) + (16008 * std::pow(N, 5) - 64052 * std::pow(N, 4) + 96094 * std::pow(N, 3) - 68071 * std::pow(N, 2) + 23024 * N - 3003) * std::pow(m, 11) - 2 * (15976 * std::pow(N, 5) - 64004 * std::pow(N, 4) + 96086 * std::pow(N, 3) - 68083 * std::pow(N, 2) + 23028 * N - 3003) * std::pow(m, 9) + (47696 * std::pow(N, 5) - 191668 * std::pow(N, 4) + 288208 * std::pow(N, 3) - 204341 * std::pow(N, 2) + 69114 * N - 9009) * std::pow(m, 8) - 8 * (6752 * std::pow(N, 5) - 27292 * std::pow(N, 4) + 41164 * std::pow(N, 3) - 29219 * std::pow(N, 2) + 9882 * N - 1287) * std::pow(m, 7) + 7 * (6656 * std::pow(N, 5) - 27160 * std::pow(N, 4) + 41164 * std::pow(N, 3) - 29270 * std::pow(N, 2) + 9897 * N - 1287) * std::pow(m, 6) - 14 * (2176 * std::pow(N, 5) - 9000 * std::pow(N, 4) + 13732 * std::pow(N, 3) - 9786 * std::pow(N, 2) + 3307 * N - 429) * std::pow(m, 5) + 16 * std::pow(N, 5) + 7 * (2128 * std::pow(N, 5) - 8952 * std::pow(N, 4) + 13768 * std::pow(N, 3) - 9834 * std::pow(N, 2) + 3319 * N - 429) * std::pow(m, 4) - 68 * std::pow(N, 4) - 4 * (1328 * std::pow(N, 5) - 5688 * std::pow(N, 4) + 8816 * std::pow(N, 3) - 6306 * std::pow(N, 2) + 2123 * N - 273) * std::pow(m, 3) + 104 * std::pow(N, 3) + (1320 * std::pow(N, 5) - 5732 * std::pow(N, 4) + 8926 * std::pow(N, 3) - 6379 * std::pow(N, 2) + 2138 * N - 273) * std::pow(m, 2) - 73 * std::pow(N, 2) - 2 * (104 * std::pow(N, 5) - 452 * std::pow(N, 4) + 702 * std::pow(N, 3) - 499 * std::pow(N, 2) + 166 * N - 21) * m + 24 * N - 3) * std::pow(r, 3) + 36 * std::pow(N, 3) + (320 * std::pow(N, 5) - 3232 * std::pow(N, 4) + 6410 * std::pow(N, 3) - 5255 * std::pow(N, 2) + 1958 * N - 273) * std::pow(m, 2) + (3 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 14) - 42 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 13) + 273 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 12) - 1092 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 11) + (48016 * std::pow(N, 5) - 192148 * std::pow(N, 4) + 288272 * std::pow(N, 3) - 204205 * std::pow(N, 2) + 69070 * N - 9009) * std::pow(m, 10) - 2 * (47888 * std::pow(N, 5) - 191972 * std::pow(N, 4) + 288208 * std::pow(N, 3) - 204209 * std::pow(N, 2) + 69074 * N - 9009) * std::pow(m, 9) + (142720 * std::pow(N, 5) - 574624 * std::pow(N, 4) + 864158 * std::pow(N, 3) - 612659 * std::pow(N, 2) + 207252 * N - 27027) * std::pow(m, 8) - 8 * (20128 * std::pow(N, 5) - 81736 * std::pow(N, 4) + 123326 * std::pow(N, 3) - 87533 * std::pow(N, 2) + 29616 * N - 3861) * std::pow(m, 7) + 7 * (19696 * std::pow(N, 5) - 81160 * std::pow(N, 4) + 123128 * std::pow(N, 3) - 87554 * std::pow(N, 2) + 29631 * N - 3861) * std::pow(m, 6) - 14 * (6352 * std::pow(N, 5) - 26776 * std::pow(N, 4) + 40952 * std::pow(N, 3) - 29198 * std::pow(N, 2) + 9885 * N - 1287) * std::pow(m, 5) + 32 * std::pow(N, 5) + 7 * (6064 * std::pow(N, 5) - 26416 * std::pow(N, 4) + 40844 * std::pow(N, 3) - 29222 * std::pow(N, 2) + 9897 * N - 1287) * std::pow(m, 4) - 176 * std::pow(N, 4) - 4 * (3632 * std::pow(N, 5) - 16544 * std::pow(N, 4) + 25924 * std::pow(N, 3) - 18622 * std::pow(N, 2) + 6309 * N - 819) * std::pow(m, 3) + 286 * std::pow(N, 3) + (3376 * std::pow(N, 5) - 16276 * std::pow(N, 4) + 25880 * std::pow(N, 3) - 18665 * std::pow(N, 2) + 6324 * N - 819) * std::pow(m, 2) - 207 * std::pow(N, 2) - 2 * (240 * std::pow(N, 5) - 1236 * std::pow(N, 4) + 1992 * std::pow(N, 3) - 1441 * std::pow(N, 2) + 488 * N - 63) * m + 70 * N - 9) * std::pow(r, 2) - 45 * std::pow(N, 2) + 2 * (144 * std::pow(N, 4) - 394 * std::pow(N, 3) + 367 * std::pow(N, 2) - 146 * N + 21) * m - (3 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 14) - 42 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 13) + 273 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 12) - 1092 * (16 * std::pow(N, 5) - 64 * std::pow(N, 4) + 96 * std::pow(N, 3) - 68 * std::pow(N, 2) + 23 * N - 3) * std::pow(m, 11) + (48000 * std::pow(N, 5) - 192116 * std::pow(N, 4) + 288236 * std::pow(N, 3) - 204185 * std::pow(N, 2) + 69066 * N - 9009) * std::pow(m, 10) - 2 * (47808 * std::pow(N, 5) - 191812 * std::pow(N, 4) + 288028 * std::pow(N, 3) - 204109 * std::pow(N, 2) + 69054 * N - 9009) * std::pow(m, 9) + 3 * (47328 * std::pow(N, 5) - 191048 * std::pow(N, 4) + 287502 * std::pow(N, 3) - 203917 * std::pow(N, 2) + 69024 * N - 9009) * std::pow(m, 8) - 24 * (6624 * std::pow(N, 5) - 27072 * std::pow(N, 4) + 40918 * std::pow(N, 3) - 29075 * std::pow(N, 2) + 9852 * N - 1287) * std::pow(m, 7) + 21 * (6384 * std::pow(N, 5) - 26680 * std::pow(N, 4) + 40640 * std::pow(N, 3) - 28974 * std::pow(N, 2) + 9837 * N - 1287) * std::pow(m, 6) - 42 * (2000 * std::pow(N, 5) - 8680 * std::pow(N, 4) + 13392 * std::pow(N, 3) - 9602 * std::pow(N, 2) + 3271 * N - 429) * std::pow(m, 5) + 21 * (1808 * std::pow(N, 5) - 8352 * std::pow(N, 4) + 13148 * std::pow(N, 3) - 9514 * std::pow(N, 2) + 3259 * N - 429) * std::pow(m, 4) - 104 * std::pow(N, 4) - 12 * (976 * std::pow(N, 5) - 5008 * std::pow(N, 4) + 8132 * std::pow(N, 3) - 5970 * std::pow(N, 2) + 2063 * N - 273) * std::pow(m, 3) + 218 * std::pow(N, 3) + 3 * (736 * std::pow(N, 5) - 4572 * std::pow(N, 4) + 7788 * std::pow(N, 3) - 5847 * std::pow(N, 2) + 2048 * N - 273) * std::pow(m, 2) - 179 * std::pow(N, 2) - 2 * (96 * std::pow(N, 5) - 916 * std::pow(N, 4) + 1684 * std::pow(N, 3) - 1309 * std::pow(N, 2) + 468 * N - 63) * m + 66 * N - 9) * r + 20 * N - 3) / ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29824 * std::pow(N, 7) - 165248 * std::pow(N, 6) + 378696 * std::pow(N, 5) - 466348 * std::pow(N, 4) + 334102 * std::pow(N, 3) - 139653 * std::pow(N, 2) + 31630 * N - 3003) * std::pow(m, 10) - 2 * (28032 * std::pow(N, 7) - 156736 * std::pow(N, 6) + 362280 * std::pow(N, 5) - 449948 * std::pow(N, 4) + 325166 * std::pow(N, 3) - 137133 * std::pow(N, 2) + 31342 * N - 3003) * std::pow(m, 9) + (76032 * std::pow(N, 7) - 433088 * std::pow(N, 6) + 1016944 * std::pow(N, 5) - 1281244 * std::pow(N, 4) + 938600 * std::pow(N, 3) - 401091 * std::pow(N, 2) + 92856 * N - 9009) * std::pow(m, 8) - 8 * (9216 * std::pow(N, 7) - 54464 * std::pow(N, 6) + 131584 * std::pow(N, 5) - 169780 * std::pow(N, 4) + 127028 * std::pow(N, 3) - 55341 * std::pow(N, 2) + 13044 * N - 1287) * std::pow(m, 7) + (49152 * std::pow(N, 7) - 312576 * std::pow(N, 6) + 794936 * std::pow(N, 5) - 1066988 * std::pow(N, 4) + 824738 * std::pow(N, 3) - 369537 * std::pow(N, 2) + 89292 * N - 9009) * std::pow(m, 6) - 2 * (10240 * std::pow(N, 7) - 76032 * std::pow(N, 6) + 212520 * std::pow(N, 5) - 304676 * std::pow(N, 4) + 247654 * std::pow(N, 3) - 115579 * std::pow(N, 2) + 28900 * N - 3003) * std::pow(m, 5) + (4096 * std::pow(N, 7) - 45056 * std::pow(N, 6) + 151696 * std::pow(N, 5) - 243524 * std::pow(N, 4) + 214156 * std::pow(N, 3) - 106063 * std::pow(N, 2) + 27802 * N - 3003) * std::pow(m, 4) + 4 * (1536 * std::pow(N, 6) - 8032 * std::pow(N, 5) + 15884 * std::pow(N, 4) - 15864 * std::pow(N, 3) + 8581 * std::pow(N, 2) - 2402 * N + 273) * std::pow(m, 3) - ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29888 * std::pow(N, 7) - 165472 * std::pow(N, 6) + 379032 * std::pow(N, 5) - 466628 * std::pow(N, 4) + 334238 * std::pow(N, 3) - 139689 * std::pow(N, 2) + 31634 * N - 3003) * std::pow(m, 10) - 2 * (28352 * std::pow(N, 7) - 157856 * std::pow(N, 6) + 363960 * std::pow(N, 5) - 451348 * std::pow(N, 4) + 325846 * std::pow(N, 3) - 137313 * std::pow(N, 2) + 31362 * N - 3003) * std::pow(m, 9) + (78848 * std::pow(N, 7) - 442976 * std::pow(N, 6) + 1031824 * std::pow(N, 5) - 1293684 * std::pow(N, 4) + 944664 * std::pow(N, 3) - 402703 * std::pow(N, 2) + 93036 * N - 9009) * std::pow(m, 8) - 8 * (10112 * std::pow(N, 7) - 57632 * std::pow(N, 6) + 136384 * std::pow(N, 5) - 173820 * std::pow(N, 4) + 129012 * std::pow(N, 3) - 55873 * std::pow(N, 2) + 13104 * N - 1287) * std::pow(m, 7) + 7 * (8672 * std::pow(N, 7) - 50576 * std::pow(N, 6) + 122656 * std::pow(N, 5) - 160176 * std::pow(N, 4) + 121674 * std::pow(N, 3) - 53839 * std::pow(N, 2) + 12876 * N - 1287) * std::pow(m, 6) - 2 * (16224 * std::pow(N, 7) - 98256 * std::pow(N, 6) + 247584 * std::pow(N, 5) - 335248 * std::pow(N, 4) + 263186 * std::pow(N, 3) - 119891 * std::pow(N, 2) + 29404 * N - 3003) * std::pow(m, 5) + 16 * std::pow(N, 5) + (11712 * std::pow(N, 7) - 75696 * std::pow(N, 6) + 202800 * std::pow(N, 5) - 290008 * std::pow(N, 4) + 238616 * std::pow(N, 3) - 113063 * std::pow(N, 2) + 28642 * N - 3003) * std::pow(m, 4) - 68 * std::pow(N, 4) - 4 * (640 * std::pow(N, 7) - 4688 * std::pow(N, 6) + 13936 * std::pow(N, 5) - 21688 * std::pow(N, 4) + 19100 * std::pow(N, 3) - 9549 * std::pow(N, 2) + 2522 * N - 273) * std::pow(m, 3) + 104 * std::pow(N, 3) + (256 * std::pow(N, 7) - 2560 * std::pow(N, 6) + 9272 * std::pow(N, 5) - 16580 * std::pow(N, 4) + 16162 * std::pow(N, 3) - 8703 * std::pow(N, 2) + 2426 * N - 273) * std::pow(m, 2) - 73 * std::pow(N, 2) + 2 * (64 * std::pow(N, 6) - 376 * std::pow(N, 5) + 868 * std::pow(N, 4) - 994 * std::pow(N, 3) + 595 * std::pow(N, 2) - 178 * N + 21) * m + 24 * N - 3) * std::pow(r, 3) + 36 * std::pow(N, 3) + (3008 * std::pow(N, 5) - 9456 * std::pow(N, 4) + 11790 * std::pow(N, 3) - 7307 * std::pow(N, 2) + 2246 * N - 273) * std::pow(m, 2) + (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89632 * std::pow(N, 7) - 496320 * std::pow(N, 6) + 1136968 * std::pow(N, 5) - 1399780 * std::pow(N, 4) + 1002660 * std::pow(N, 3) - 419051 * std::pow(N, 2) + 94900 * N - 9009) * std::pow(m, 10) - 2 * (84896 * std::pow(N, 7) - 473088 * std::pow(N, 6) + 1091240 * std::pow(N, 5) - 1353524 * std::pow(N, 4) + 977268 * std::pow(N, 3) - 411859 * std::pow(N, 2) + 94076 * N - 9009) * std::pow(m, 9) + (235136 * std::pow(N, 7) - 1324688 * std::pow(N, 6) + 3089800 * std::pow(N, 5) - 3876432 * std::pow(N, 4) + 2831586 * std::pow(N, 3) - 1207393 * std::pow(N, 2) + 279018 * N - 27027) * std::pow(m, 8) - 8 * (29888 * std::pow(N, 7) - 171536 * std::pow(N, 6) + 407320 * std::pow(N, 5) - 519960 * std::pow(N, 4) + 386250 * std::pow(N, 3) - 167383 * std::pow(N, 2) + 39282 * N - 3861) * std::pow(m, 7) + (176320 * std::pow(N, 7) - 1044240 * std::pow(N, 6) + 2551416 * std::pow(N, 5) - 3343560 * std::pow(N, 4) + 2544490 * std::pow(N, 3) - 1127371 * std::pow(N, 2) + 269976 * N - 27027) * std::pow(m, 6) - 2 * (45632 * std::pow(N, 7) - 285104 * std::pow(N, 6) + 729256 * std::pow(N, 5) - 994392 * std::pow(N, 4) + 783438 * std::pow(N, 3) - 357769 * std::pow(N, 2) + 87960 * N - 9009) * std::pow(m, 5) + 32 * std::pow(N, 5) + (31104 * std::pow(N, 7) - 213440 * std::pow(N, 6) + 588488 * std::pow(N, 5) - 852736 * std::pow(N, 4) + 706248 * std::pow(N, 3) - 336109 * std::pow(N, 2) + 85506 * N - 9009) * std::pow(m, 4) - 176 * std::pow(N, 4) - 4 * (1536 * std::pow(N, 7) - 12576 * std::pow(N, 6) + 39448 * std::pow(N, 5) - 62896 * std::pow(N, 4) + 56036 * std::pow(N, 3) - 28223 * std::pow(N, 2) + 7506 * N - 819) * std::pow(m, 3) + 286 * std::pow(N, 3) + (512 * std::pow(N, 7) - 6272 * std::pow(N, 6) + 25168 * std::pow(N, 5) - 47052 * std::pow(N, 4) + 46788 * std::pow(N, 3) - 25501 * std::pow(N, 2) + 7188 * N - 819) * std::pow(m, 2) - 207 * std::pow(N, 2) + 2 * (128 * std::pow(N, 6) - 944 * std::pow(N, 5) + 2380 * std::pow(N, 4) - 2820 * std::pow(N, 3) + 1721 * std::pow(N, 2) - 524 * N + 63) * m + 70 * N - 9) * std::pow(r, 2) - 45 * std::pow(N, 2) + 2 * (288 * std::pow(N, 4) - 574 * std::pow(N, 3) + 447 * std::pow(N, 2) - 158 * N + 21) * m - (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89568 * std::pow(N, 7) - 496096 * std::pow(N, 6) + 1136632 * std::pow(N, 5) - 1399500 * std::pow(N, 4) + 1002524 * std::pow(N, 3) - 419015 * std::pow(N, 2) + 94896 * N - 9009) * std::pow(m, 10) - 2 * (84576 * std::pow(N, 7) - 471968 * std::pow(N, 6) + 1089560 * std::pow(N, 5) - 1352124 * std::pow(N, 4) + 976588 * std::pow(N, 3) - 411679 * std::pow(N, 2) + 94056 * N - 9009) * std::pow(m, 9) + (232320 * std::pow(N, 7) - 1314800 * std::pow(N, 6) + 3074920 * std::pow(N, 5) - 3863992 * std::pow(N, 4) + 2825522 * std::pow(N, 3) - 1205781 * std::pow(N, 2) + 278838 * N - 27027) * std::pow(m, 8) - 8 * (28992 * std::pow(N, 7) - 168368 * std::pow(N, 6) + 402520 * std::pow(N, 5) - 515920 * std::pow(N, 4) + 384266 * std::pow(N, 3) - 166851 * std::pow(N, 2) + 39222 * N - 3861) * std::pow(m, 7) + (164736 * std::pow(N, 7) - 1002736 * std::pow(N, 6) + 2487736 * std::pow(N, 5) - 3289312 * std::pow(N, 4) + 2517510 * std::pow(N, 3) - 1120035 * std::pow(N, 2) + 269136 * N - 27027) * std::pow(m, 6) - 2 * (39552 * std::pow(N, 7) - 262736 * std::pow(N, 6) + 694120 * std::pow(N, 5) - 963808 * std::pow(N, 4) + 767906 * std::pow(N, 3) - 353457 * std::pow(N, 2) + 87456 * N - 9009) * std::pow(m, 5) + (23040 * std::pow(N, 7) - 182112 * std::pow(N, 6) + 537032 * std::pow(N, 5) - 806192 * std::pow(N, 4) + 681788 * std::pow(N, 3) - 329109 * std::pow(N, 2) + 84666 * N - 9009) * std::pow(m, 4) - 104 * std::pow(N, 4) - 12 * (256 * std::pow(N, 7) - 3072 * std::pow(N, 6) + 11144 * std::pow(N, 5) - 19024 * std::pow(N, 4) + 17600 * std::pow(N, 3) - 9085 * std::pow(N, 2) + 2462 * N - 273) * std::pow(m, 3) + 218 * std::pow(N, 3) - (3200 * std::pow(N, 6) - 18592 * std::pow(N, 5) + 39868 * std::pow(N, 4) - 42416 * std::pow(N, 3) + 24105 * std::pow(N, 2) - 7008 * N + 819) * std::pow(m, 2) - 179 * std::pow(N, 2) - 2 * (512 * std::pow(N, 5) - 1788 * std::pow(N, 4) + 2400 * std::pow(N, 3) - 1573 * std::pow(N, 2) + 504 * N - 63) * m + 66 * N - 9) * r + 20 * N - 3);
            return result;
        }

        static double calc_deltaAB(double N, double m, double r)
        {
            double result = ((8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 14) - 14 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 13) + (736 * std::pow(N, 5) - 4028 * std::pow(N, 4) + 7124 * std::pow(N, 3) - 5563 * std::pow(N, 2) + 2004 * N - 273) * std::pow(m, 12) - 4 * (752 * std::pow(N, 5) - 4076 * std::pow(N, 4) + 7176 * std::pow(N, 3) - 5587 * std::pow(N, 2) + 2008 * N - 273) * std::pow(m, 11) + (8520 * std::pow(N, 5) - 45596 * std::pow(N, 4) + 79758 * std::pow(N, 3) - 61833 * std::pow(N, 2) + 22150 * N - 3003) * std::pow(m, 10) - 2 * (8808 * std::pow(N, 5) - 46524 * std::pow(N, 4) + 80758 * std::pow(N, 3) - 62281 * std::pow(N, 2) + 22222 * N - 3003) * std::pow(m, 9) + (27248 * std::pow(N, 5) - 142532 * std::pow(N, 4) + 245452 * std::pow(N, 3) - 188215 * std::pow(N, 2) + 66876 * N - 9009) * std::pow(m, 8) - 8 * (3968 * std::pow(N, 5) - 20732 * std::pow(N, 4) + 35464 * std::pow(N, 3) - 27049 * std::pow(N, 2) + 9576 * N - 1287) * std::pow(m, 7) + (27616 * std::pow(N, 5) - 146472 * std::pow(N, 4) + 249802 * std::pow(N, 3) - 189847 * std::pow(N, 2) + 67074 * N - 9009) * std::pow(m, 6) - 2 * (8736 * std::pow(N, 5) - 48408 * std::pow(N, 4) + 82942 * std::pow(N, 3) - 63069 * std::pow(N, 2) + 22310 * N - 3003) * std::pow(m, 5) + (7616 * std::pow(N, 5) - 46464 * std::pow(N, 4) + 81148 * std::pow(N, 3) - 62241 * std::pow(N, 2) + 22172 * N - 3003) * std::pow(m, 4) - 4 * (512 * std::pow(N, 5) - 3824 * std::pow(N, 4) + 7000 * std::pow(N, 3) - 5499 * std::pow(N, 2) + 1992 * N - 273) * std::pow(m, 3) - ((8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 14) - 14 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 13) + (736 * std::pow(N, 5) - 4028 * std::pow(N, 4) + 7124 * std::pow(N, 3) - 5563 * std::pow(N, 2) + 2004 * N - 273) * std::pow(m, 12) - 4 * (752 * std::pow(N, 5) - 4076 * std::pow(N, 4) + 7176 * std::pow(N, 3) - 5587 * std::pow(N, 2) + 2008 * N - 273) * std::pow(m, 11) + 11 * (776 * std::pow(N, 5) - 4148 * std::pow(N, 4) + 7254 * std::pow(N, 3) - 5623 * std::pow(N, 2) + 2014 * N - 273) * std::pow(m, 10) - 22 * (808 * std::pow(N, 5) - 4244 * std::pow(N, 4) + 7358 * std::pow(N, 3) - 5671 * std::pow(N, 2) + 2022 * N - 273) * std::pow(m, 9) + 33 * (848 * std::pow(N, 5) - 4364 * std::pow(N, 4) + 7488 * std::pow(N, 3) - 5731 * std::pow(N, 2) + 2032 * N - 273) * std::pow(m, 8) - 264 * (128 * std::pow(N, 5) - 644 * std::pow(N, 4) + 1092 * std::pow(N, 3) - 829 * std::pow(N, 2) + 292 * N - 39) * std::pow(m, 7) + 231 * (136 * std::pow(N, 5) - 668 * std::pow(N, 4) + 1118 * std::pow(N, 3) - 841 * std::pow(N, 2) + 294 * N - 39) * std::pow(m, 6) - 22 * (1016 * std::pow(N, 5) - 4868 * std::pow(N, 4) + 8034 * std::pow(N, 3) - 5983 * std::pow(N, 2) + 2074 * N - 273) * std::pow(m, 5) + 16 * std::pow(N, 5) + 11 * (1088 * std::pow(N, 5) - 5084 * std::pow(N, 4) + 8268 * std::pow(N, 3) - 6091 * std::pow(N, 2) + 2092 * N - 273) * std::pow(m, 4) - 68 * std::pow(N, 4) - 4 * (1168 * std::pow(N, 5) - 5324 * std::pow(N, 4) + 8528 * std::pow(N, 3) - 6211 * std::pow(N, 2) + 2112 * N - 273) * std::pow(m, 3) + 104 * std::pow(N, 3) + (1256 * std::pow(N, 5) - 5588 * std::pow(N, 4) + 8814 * std::pow(N, 3) - 6343 * std::pow(N, 2) + 2134 * N - 273) * std::pow(m, 2) - 73 * std::pow(N, 2) - 2 * (104 * std::pow(N, 5) - 452 * std::pow(N, 4) + 702 * std::pow(N, 3) - 499 * std::pow(N, 2) + 166 * N - 21) * m + 24 * N - 3) * std::pow(r, 3) + 36 * std::pow(N, 3) + (256 * std::pow(N, 5) - 3088 * std::pow(N, 4) + 6298 * std::pow(N, 3) - 5219 * std::pow(N, 2) + 1954 * N - 273) * std::pow(m, 2) + (3 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 14) - 42 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 13) + 3 * (736 * std::pow(N, 5) - 4028 * std::pow(N, 4) + 7124 * std::pow(N, 3) - 5563 * std::pow(N, 2) + 2004 * N - 273) * std::pow(m, 12) - 12 * (752 * std::pow(N, 5) - 4076 * std::pow(N, 4) + 7176 * std::pow(N, 3) - 5587 * std::pow(N, 2) + 2008 * N - 273) * std::pow(m, 11) + (25600 * std::pow(N, 5) - 136876 * std::pow(N, 4) + 239372 * std::pow(N, 3) - 185551 * std::pow(N, 2) + 66460 * N - 9009) * std::pow(m, 10) - 2 * (26624 * std::pow(N, 5) - 140012 * std::pow(N, 4) + 242764 * std::pow(N, 3) - 187103 * std::pow(N, 2) + 66716 * N - 9009) * std::pow(m, 9) + (83584 * std::pow(N, 5) - 431656 * std::pow(N, 4) + 740846 * std::pow(N, 3) - 567005 * std::pow(N, 2) + 201078 * N - 27027) * std::pow(m, 8) - 8 * (12544 * std::pow(N, 5) - 63616 * std::pow(N, 4) + 107942 * std::pow(N, 3) - 81947 * std::pow(N, 2) + 28878 * N - 3861) * std::pow(m, 7) + 7 * (13192 * std::pow(N, 5) - 65812 * std::pow(N, 4) + 110318 * std::pow(N, 3) - 83003 * std::pow(N, 2) + 29046 * N - 3861) * std::pow(m, 6) - 2 * (32296 * std::pow(N, 5) - 159076 * std::pow(N, 4) + 263414 * std::pow(N, 3) - 196319 * std::pow(N, 2) + 68190 * N - 9009) * std::pow(m, 5) + 32 * std::pow(N, 5) + (33664 * std::pow(N, 5) - 164692 * std::pow(N, 4) + 269624 * std::pow(N, 3) - 199043 * std::pow(N, 2) + 68616 * N - 9009) * std::pow(m, 4) - 176 * std::pow(N, 4) - 4 * (3152 * std::pow(N, 5) - 15452 * std::pow(N, 4) + 25060 * std::pow(N, 3) - 18337 * std::pow(N, 2) + 6276 * N - 819) * std::pow(m, 3) + 286 * std::pow(N, 3) + (3184 * std::pow(N, 5) - 15844 * std::pow(N, 4) + 25544 * std::pow(N, 3) - 18557 * std::pow(N, 2) + 6312 * N - 819) * std::pow(m, 2) - 207 * std::pow(N, 2) - 2 * (240 * std::pow(N, 5) - 1236 * std::pow(N, 4) + 1992 * std::pow(N, 3) - 1441 * std::pow(N, 2) + 488 * N - 63) * m + 70 * N - 9) * std::pow(r, 2) - 45 * std::pow(N, 2) + 2 * (144 * std::pow(N, 4) - 394 * std::pow(N, 3) + 367 * std::pow(N, 2) - 146 * N + 21) * m - (3 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 14) - 42 * (8 * std::pow(N, 5) - 44 * std::pow(N, 4) + 78 * std::pow(N, 3) - 61 * std::pow(N, 2) + 22 * N - 3) * std::pow(m, 13) + 3 * (736 * std::pow(N, 5) - 4028 * std::pow(N, 4) + 7124 * std::pow(N, 3) - 5563 * std::pow(N, 2) + 2004 * N - 273) * std::pow(m, 12) - 12 * (752 * std::pow(N, 5) - 4076 * std::pow(N, 4) + 7176 * std::pow(N, 3) - 5587 * std::pow(N, 2) + 2008 * N - 273) * std::pow(m, 11) + (25584 * std::pow(N, 5) - 136844 * std::pow(N, 4) + 239336 * std::pow(N, 3) - 185531 * std::pow(N, 2) + 66456 * N - 9009) * std::pow(m, 10) - 2 * (26544 * std::pow(N, 5) - 139852 * std::pow(N, 4) + 242584 * std::pow(N, 3) - 187003 * std::pow(N, 2) + 66696 * N - 9009) * std::pow(m, 9) + 3 * (27616 * std::pow(N, 5) - 143392 * std::pow(N, 4) + 246398 * std::pow(N, 3) - 188699 * std::pow(N, 2) + 66966 * N - 9009) * std::pow(m, 8) - 24 * (4096 * std::pow(N, 5) - 21032 * std::pow(N, 4) + 35790 * std::pow(N, 3) - 27213 * std::pow(N, 2) + 9606 * N - 1287) * std::pow(m, 7) + 21 * (4216 * std::pow(N, 5) - 21564 * std::pow(N, 4) + 36370 * std::pow(N, 3) - 27457 * std::pow(N, 2) + 9642 * N - 1287) * std::pow(m, 6) - 6 * (9944 * std::pow(N, 5) - 51308 * std::pow(N, 4) + 85994 * std::pow(N, 3) - 64525 * std::pow(N, 2) + 22562 * N - 3003) * std::pow(m, 5) + 3 * (9728 * std::pow(N, 5) - 51724 * std::pow(N, 4) + 86608 * std::pow(N, 3) - 64761 * std::pow(N, 2) + 22592 * N - 3003) * std::pow(m, 4) - 104 * std::pow(N, 4) - 12 * (816 * std::pow(N, 5) - 4644 * std::pow(N, 4) + 7844 * std::pow(N, 3) - 5875 * std::pow(N, 2) + 2052 * N - 273) * std::pow(m, 3) + 218 * std::pow(N, 3) + 3 * (672 * std::pow(N, 5) - 4428 * std::pow(N, 4) + 7676 * std::pow(N, 3) - 5811 * std::pow(N, 2) + 2044 * N - 273) * std::pow(m, 2) - 179 * std::pow(N, 2) - 2 * (96 * std::pow(N, 5) - 916 * std::pow(N, 4) + 1684 * std::pow(N, 3) - 1309 * std::pow(N, 2) + 468 * N - 63) * m + 66 * N - 9) * r + 20 * N - 3) / ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29824 * std::pow(N, 7) - 165248 * std::pow(N, 6) + 378696 * std::pow(N, 5) - 466348 * std::pow(N, 4) + 334102 * std::pow(N, 3) - 139653 * std::pow(N, 2) + 31630 * N - 3003) * std::pow(m, 10) - 2 * (28032 * std::pow(N, 7) - 156736 * std::pow(N, 6) + 362280 * std::pow(N, 5) - 449948 * std::pow(N, 4) + 325166 * std::pow(N, 3) - 137133 * std::pow(N, 2) + 31342 * N - 3003) * std::pow(m, 9) + (76032 * std::pow(N, 7) - 433088 * std::pow(N, 6) + 1016944 * std::pow(N, 5) - 1281244 * std::pow(N, 4) + 938600 * std::pow(N, 3) - 401091 * std::pow(N, 2) + 92856 * N - 9009) * std::pow(m, 8) - 8 * (9216 * std::pow(N, 7) - 54464 * std::pow(N, 6) + 131584 * std::pow(N, 5) - 169780 * std::pow(N, 4) + 127028 * std::pow(N, 3) - 55341 * std::pow(N, 2) + 13044 * N - 1287) * std::pow(m, 7) + (49152 * std::pow(N, 7) - 312576 * std::pow(N, 6) + 794936 * std::pow(N, 5) - 1066988 * std::pow(N, 4) + 824738 * std::pow(N, 3) - 369537 * std::pow(N, 2) + 89292 * N - 9009) * std::pow(m, 6) - 2 * (10240 * std::pow(N, 7) - 76032 * std::pow(N, 6) + 212520 * std::pow(N, 5) - 304676 * std::pow(N, 4) + 247654 * std::pow(N, 3) - 115579 * std::pow(N, 2) + 28900 * N - 3003) * std::pow(m, 5) + (4096 * std::pow(N, 7) - 45056 * std::pow(N, 6) + 151696 * std::pow(N, 5) - 243524 * std::pow(N, 4) + 214156 * std::pow(N, 3) - 106063 * std::pow(N, 2) + 27802 * N - 3003) * std::pow(m, 4) + 4 * (1536 * std::pow(N, 6) - 8032 * std::pow(N, 5) + 15884 * std::pow(N, 4) - 15864 * std::pow(N, 3) + 8581 * std::pow(N, 2) - 2402 * N + 273) * std::pow(m, 3) - ((32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 14 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 4 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (29888 * std::pow(N, 7) - 165472 * std::pow(N, 6) + 379032 * std::pow(N, 5) - 466628 * std::pow(N, 4) + 334238 * std::pow(N, 3) - 139689 * std::pow(N, 2) + 31634 * N - 3003) * std::pow(m, 10) - 2 * (28352 * std::pow(N, 7) - 157856 * std::pow(N, 6) + 363960 * std::pow(N, 5) - 451348 * std::pow(N, 4) + 325846 * std::pow(N, 3) - 137313 * std::pow(N, 2) + 31362 * N - 3003) * std::pow(m, 9) + (78848 * std::pow(N, 7) - 442976 * std::pow(N, 6) + 1031824 * std::pow(N, 5) - 1293684 * std::pow(N, 4) + 944664 * std::pow(N, 3) - 402703 * std::pow(N, 2) + 93036 * N - 9009) * std::pow(m, 8) - 8 * (10112 * std::pow(N, 7) - 57632 * std::pow(N, 6) + 136384 * std::pow(N, 5) - 173820 * std::pow(N, 4) + 129012 * std::pow(N, 3) - 55873 * std::pow(N, 2) + 13104 * N - 1287) * std::pow(m, 7) + 7 * (8672 * std::pow(N, 7) - 50576 * std::pow(N, 6) + 122656 * std::pow(N, 5) - 160176 * std::pow(N, 4) + 121674 * std::pow(N, 3) - 53839 * std::pow(N, 2) + 12876 * N - 1287) * std::pow(m, 6) - 2 * (16224 * std::pow(N, 7) - 98256 * std::pow(N, 6) + 247584 * std::pow(N, 5) - 335248 * std::pow(N, 4) + 263186 * std::pow(N, 3) - 119891 * std::pow(N, 2) + 29404 * N - 3003) * std::pow(m, 5) + 16 * std::pow(N, 5) + (11712 * std::pow(N, 7) - 75696 * std::pow(N, 6) + 202800 * std::pow(N, 5) - 290008 * std::pow(N, 4) + 238616 * std::pow(N, 3) - 113063 * std::pow(N, 2) + 28642 * N - 3003) * std::pow(m, 4) - 68 * std::pow(N, 4) - 4 * (640 * std::pow(N, 7) - 4688 * std::pow(N, 6) + 13936 * std::pow(N, 5) - 21688 * std::pow(N, 4) + 19100 * std::pow(N, 3) - 9549 * std::pow(N, 2) + 2522 * N - 273) * std::pow(m, 3) + 104 * std::pow(N, 3) + (256 * std::pow(N, 7) - 2560 * std::pow(N, 6) + 9272 * std::pow(N, 5) - 16580 * std::pow(N, 4) + 16162 * std::pow(N, 3) - 8703 * std::pow(N, 2) + 2426 * N - 273) * std::pow(m, 2) - 73 * std::pow(N, 2) + 2 * (64 * std::pow(N, 6) - 376 * std::pow(N, 5) + 868 * std::pow(N, 4) - 994 * std::pow(N, 3) + 595 * std::pow(N, 2) - 178 * N + 21) * m + 24 * N - 3) * std::pow(r, 3) + 36 * std::pow(N, 3) + (3008 * std::pow(N, 5) - 9456 * std::pow(N, 4) + 11790 * std::pow(N, 3) - 7307 * std::pow(N, 2) + 2246 * N - 273) * std::pow(m, 2) + (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89632 * std::pow(N, 7) - 496320 * std::pow(N, 6) + 1136968 * std::pow(N, 5) - 1399780 * std::pow(N, 4) + 1002660 * std::pow(N, 3) - 419051 * std::pow(N, 2) + 94900 * N - 9009) * std::pow(m, 10) - 2 * (84896 * std::pow(N, 7) - 473088 * std::pow(N, 6) + 1091240 * std::pow(N, 5) - 1353524 * std::pow(N, 4) + 977268 * std::pow(N, 3) - 411859 * std::pow(N, 2) + 94076 * N - 9009) * std::pow(m, 9) + (235136 * std::pow(N, 7) - 1324688 * std::pow(N, 6) + 3089800 * std::pow(N, 5) - 3876432 * std::pow(N, 4) + 2831586 * std::pow(N, 3) - 1207393 * std::pow(N, 2) + 279018 * N - 27027) * std::pow(m, 8) - 8 * (29888 * std::pow(N, 7) - 171536 * std::pow(N, 6) + 407320 * std::pow(N, 5) - 519960 * std::pow(N, 4) + 386250 * std::pow(N, 3) - 167383 * std::pow(N, 2) + 39282 * N - 3861) * std::pow(m, 7) + (176320 * std::pow(N, 7) - 1044240 * std::pow(N, 6) + 2551416 * std::pow(N, 5) - 3343560 * std::pow(N, 4) + 2544490 * std::pow(N, 3) - 1127371 * std::pow(N, 2) + 269976 * N - 27027) * std::pow(m, 6) - 2 * (45632 * std::pow(N, 7) - 285104 * std::pow(N, 6) + 729256 * std::pow(N, 5) - 994392 * std::pow(N, 4) + 783438 * std::pow(N, 3) - 357769 * std::pow(N, 2) + 87960 * N - 9009) * std::pow(m, 5) + 32 * std::pow(N, 5) + (31104 * std::pow(N, 7) - 213440 * std::pow(N, 6) + 588488 * std::pow(N, 5) - 852736 * std::pow(N, 4) + 706248 * std::pow(N, 3) - 336109 * std::pow(N, 2) + 85506 * N - 9009) * std::pow(m, 4) - 176 * std::pow(N, 4) - 4 * (1536 * std::pow(N, 7) - 12576 * std::pow(N, 6) + 39448 * std::pow(N, 5) - 62896 * std::pow(N, 4) + 56036 * std::pow(N, 3) - 28223 * std::pow(N, 2) + 7506 * N - 819) * std::pow(m, 3) + 286 * std::pow(N, 3) + (512 * std::pow(N, 7) - 6272 * std::pow(N, 6) + 25168 * std::pow(N, 5) - 47052 * std::pow(N, 4) + 46788 * std::pow(N, 3) - 25501 * std::pow(N, 2) + 7188 * N - 819) * std::pow(m, 2) - 207 * std::pow(N, 2) + 2 * (128 * std::pow(N, 6) - 944 * std::pow(N, 5) + 2380 * std::pow(N, 4) - 2820 * std::pow(N, 3) + 1721 * std::pow(N, 2) - 524 * N + 63) * m + 70 * N - 9) * std::pow(r, 2) - 45 * std::pow(N, 2) + 2 * (288 * std::pow(N, 4) - 574 * std::pow(N, 3) + 447 * std::pow(N, 2) - 158 * N + 21) * m - (3 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 14) - 42 * (32 * std::pow(N, 7) - 176 * std::pow(N, 6) + 400 * std::pow(N, 5) - 488 * std::pow(N, 4) + 346 * std::pow(N, 3) - 143 * std::pow(N, 2) + 32 * N - 3) * std::pow(m, 13) + 3 * (2880 * std::pow(N, 7) - 15856 * std::pow(N, 6) + 36080 * std::pow(N, 5) - 44080 * std::pow(N, 4) + 31304 * std::pow(N, 3) - 12961 * std::pow(N, 2) + 2906 * N - 273) * std::pow(m, 12) - 12 * (2816 * std::pow(N, 7) - 15536 * std::pow(N, 6) + 35440 * std::pow(N, 5) - 43424 * std::pow(N, 4) + 30940 * std::pow(N, 3) - 12857 * std::pow(N, 2) + 2894 * N - 273) * std::pow(m, 11) + (89568 * std::pow(N, 7) - 496096 * std::pow(N, 6) + 1136632 * std::pow(N, 5) - 1399500 * std::pow(N, 4) + 1002524 * std::pow(N, 3) - 419015 * std::pow(N, 2) + 94896 * N - 9009) * std::pow(m, 10) - 2 * (84576 * std::pow(N, 7) - 471968 * std::pow(N, 6) + 1089560 * std::pow(N, 5) - 1352124 * std::pow(N, 4) + 976588 * std::pow(N, 3) - 411679 * std::pow(N, 2) + 94056 * N - 9009) * std::pow(m, 9) + (232320 * std::pow(N, 7) - 1314800 * std::pow(N, 6) + 3074920 * std::pow(N, 5) - 3863992 * std::pow(N, 4) + 2825522 * std::pow(N, 3) - 1205781 * std::pow(N, 2) + 278838 * N - 27027) * std::pow(m, 8) - 8 * (28992 * std::pow(N, 7) - 168368 * std::pow(N, 6) + 402520 * std::pow(N, 5) - 515920 * std::pow(N, 4) + 384266 * std::pow(N, 3) - 166851 * std::pow(N, 2) + 39222 * N - 3861) * std::pow(m, 7) + (164736 * std::pow(N, 7) - 1002736 * std::pow(N, 6) + 2487736 * std::pow(N, 5) - 3289312 * std::pow(N, 4) + 2517510 * std::pow(N, 3) - 1120035 * std::pow(N, 2) + 269136 * N - 27027) * std::pow(m, 6) - 2 * (39552 * std::pow(N, 7) - 262736 * std::pow(N, 6) + 694120 * std::pow(N, 5) - 963808 * std::pow(N, 4) + 767906 * std::pow(N, 3) - 353457 * std::pow(N, 2) + 87456 * N - 9009) * std::pow(m, 5) + (23040 * std::pow(N, 7) - 182112 * std::pow(N, 6) + 537032 * std::pow(N, 5) - 806192 * std::pow(N, 4) + 681788 * std::pow(N, 3) - 329109 * std::pow(N, 2) + 84666 * N - 9009) * std::pow(m, 4) - 104 * std::pow(N, 4) - 12 * (256 * std::pow(N, 7) - 3072 * std::pow(N, 6) + 11144 * std::pow(N, 5) - 19024 * std::pow(N, 4) + 17600 * std::pow(N, 3) - 9085 * std::pow(N, 2) + 2462 * N - 273) * std::pow(m, 3) + 218 * std::pow(N, 3) - (3200 * std::pow(N, 6) - 18592 * std::pow(N, 5) + 39868 * std::pow(N, 4) - 42416 * std::pow(N, 3) + 24105 * std::pow(N, 2) - 7008 * N + 819) * std::pow(m, 2) - 179 * std::pow(N, 2) - 2 * (512 * std::pow(N, 5) - 1788 * std::pow(N, 4) + 2400 * std::pow(N, 3) - 1573 * std::pow(N, 2) + 504 * N - 63) * m + 66 * N - 9) * r + 20 * N - 3);
            return result;
        }
    };
} // namespace unit_test

TEST_CASE_METHOD(unit_test::FunctionnalTest, "calc_phiA_test")
{
    auto result = calc_phiA(1000, 0.0001);

    REQUIRE(result == Approx(0.7142).epsilon(0.001));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "calc_phiAB_test")
{
    auto result = calc_phiAB(1000, 0.0001, 0.0001);

    REQUIRE(result == Approx(0.5490).epsilon(0.001));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "calc_gammaAB_test")
{
    auto result = calc_gammaAB(1000, 0.0001, 0.0001);

    REQUIRE(result == Approx(0.5201).epsilon(0.001));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "calc_deltaAB_test")
{
    auto result = calc_deltaAB(1000, 0.0001, 0.0001);

    REQUIRE(result == Approx(0.5160).epsilon(0.001));
}

/****************************/
/*       Algo Hudson        */
/****************************/

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_hudson_site_2_N_5000_diplo_m_0_0001_r_0_0001_rep_400000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 400000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.0001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_hudson_site_2_N_5000_diplo_m_0_00001_r_0_0001_rep_800000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 800000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.0001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.00001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_hudson_site_2_N_5000_diplo_m_0_0001_r_0_00001_rep_400000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 400000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.00001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}

/****************************/
/*      Algo Gen/Gen        */
/****************************/

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_gen/gen_one_pop_site_2_N_5000_diplo_m_0_0001_r_0_0001_rep_1200000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1200000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.0001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_gen/gen_one_pop_site_2_N_5000_diplo_m_0_00001_r_0_0001_rep_1200000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1200000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.0001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.00001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}

TEST_CASE_METHOD(unit_test::FunctionnalTest, "algo_gen/gen_one_pop_site_2_N_5000_diplo_m_0_0001_r_0_00001_rep_1600000_iam")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1600000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Prob_id_1_2_loc = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_size_per_node = std::vector<int>(1, 20);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 2;
    recomb_param.Unscaled_recomb_rate = 0.00001;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 0));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    /*****************/

    auto phiAB = calc_phiAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto gammaAB = calc_gammaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto deltaAB = calc_deltaAB(demo_param.Population_size_N, muta_param.Unscaled_mut_rate_mu, recomb_param.Unscaled_recomb_rate);
    auto PHI = (phiAB + 2 * gammaAB + deltaAB) / 4.0;

    auto m_PHI = std::get<0>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    auto v_PHI = std::get<1>(output.Prob_id_1_2_loc_res.PHI_cumul_m_v);
    REQUIRE(m_PHI == Approx(PHI).margin(0.0001));
    REQUIRE(PHI < m_PHI + 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
    REQUIRE(PHI > m_PHI - 1.96 * (std::sqrt(v_PHI / simu_param.Repetition_nbr)));
}