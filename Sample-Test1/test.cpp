#include "pch.h"
#include "..\ComputationalMechanichs\cart2dec.h"
#include "..\ComputationalMechanichs\dec2cart.h"

using namespace N;
using namespace C;

TEST(TestCaseName, TestName) {
    for (int i = 0; i < 500; i++) {
        std::vector<double> first_ = { 2e5, 0.3, 0.3, 1.5, 1.5, 0 };
        cart2dec first = cart2dec(first_[0], first_[1], first_[2], first_[3], first_[4], first_[5]);
        std::vector<double> first_trans = first.coordinates_velocities(i);
        dec2cart second = dec2cart(first_trans[0], first_trans[1], first_trans[2], first_trans[3], first_trans[4], first_trans[5]);
        std::vector<double> second_trans = second.keplerian_elements(i);
        EXPECT_NEAR(first_[0], second_trans[0], first_[0]*5e-4);
        EXPECT_NEAR(first_[1], second_trans[1], first_[1] * 5e-4);
        EXPECT_NEAR(first_[2], second_trans[2], first_[2] * 1e-4);
        EXPECT_NEAR(first_[3], second_trans[3], first_[3] * 1e-4);
        EXPECT_NEAR(first_[4], second_trans[4], first_[4] * 5e-4);
        EXPECT_NEAR(first_[5], second_trans[5], i * 0.35);
    }
}

TEST(TestCaseName1, TestName1) {
    for (int i = 0; i < 100; i++) {
        std::vector<double> first_ = { 1.5e5, 0.5, 0.2, 1, 1, 0 };
        cart2dec first = cart2dec(first_[0], first_[1], first_[2], first_[3], first_[4], first_[5]);
        std::vector<double> first_trans = first.coordinates_velocities(i);
        dec2cart second = dec2cart(first_trans[0], first_trans[1], first_trans[2], first_trans[3], first_trans[4], first_trans[5]);
        std::vector<double> second_trans = second.keplerian_elements(i);
        EXPECT_NEAR(first_[0], second_trans[0], first_[0] * 5e-4);
        EXPECT_NEAR(first_[1], second_trans[1], first_[1] * 5e-4);
        EXPECT_NEAR(first_[2], second_trans[2], first_[2] * 1e-4);
        EXPECT_NEAR(first_[3], second_trans[3], first_[3] * 1e-4);
        EXPECT_NEAR(first_[4], second_trans[4], first_[4] * 5e-4);
        EXPECT_NEAR(first_[5], second_trans[5], i * 0.35);
    }
}

TEST(TestCaseName2, TestName2) {
    for (int i = 0; i < 100; i++) {
        std::vector<double> first_ = { 1e5, 0.3, 0.2, 1, 1, 0 };
        cart2dec first = cart2dec(first_[0], first_[1], first_[2], first_[3], first_[4], first_[5]);
        std::vector<double> first_trans = first.coordinates_velocities(i);
        dec2cart second = dec2cart(first_trans[0], first_trans[1], first_trans[2], first_trans[3], first_trans[4], first_trans[5]);
        std::vector<double> second_trans = second.keplerian_elements(i);
        EXPECT_NEAR(first_[0], second_trans[0], first_[0] * 5e-4);
        EXPECT_NEAR(first_[1], second_trans[1], first_[1] * 5e-4);
        EXPECT_NEAR(first_[2], second_trans[2], first_[2] * 1e-4);
        EXPECT_NEAR(first_[3], second_trans[3], first_[3] * 1e-4);
        EXPECT_NEAR(first_[4], second_trans[4], first_[4] * 5e-4);
        EXPECT_NEAR(first_[5], second_trans[5], i * 0.35);
    }
}