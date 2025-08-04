//
// Created by PAndA on 25-07-2025.
//

#include "../include/D2Q9.hpp"
#include <iostream>
#include <array>
#include <fstream>


void d2q9_sim() {
    double alpha = 0.25;
    double omega = 1.0 / (3. * alpha + 0.5);


    std::array<std::array<double, M>, M> rho{};
    std::array<double, M> row{};
    row.fill(1.0);

    rho.fill(row);
    std::array<std::array<double, M>, M> f0, f1, f2, f3, f4, f5, f6, f7, f8;

    std::array<std::array<double, M>, M> f0_old, f1_old, f2_old, f3_old, f4_old, f5_old, f6_old, f7_old, f8_old;

    std::array<std::array<double, M>, M> feq0, feq1, feq2, feq3, feq4, feq5, feq6, feq7, feq8;

    f0.fill(row);
    f1.fill(row);
    f2.fill(row);
    f3.fill(row);
    f4.fill(row);
    f5.fill(row);
    f6.fill(row);
    f7.fill(row);
    f8.fill(row);


    feq0.fill(row);
    feq1.fill(row);
    feq2.fill(row);
    feq3.fill(row);
    feq4.fill(row);
    feq5.fill(row);
    feq6.fill(row);
    feq7.fill(row);
    feq8.fill(row);

    unsigned int ntimes = 100;
    const float tw = 1.0f;

    for (int t = 0; t < ntimes; ++t) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < M; ++i) {
                rho[i][j] = f0[i][j] + f1[i][j] + f2[i][j] + f3[i][j] + f4[i][j] + f5[i][j] + f6[i][j] + f7[i][j] + f8[
                                i][j];
            }
        }
        // ===================          1. collision            ====================
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < M; ++i) {
                feq0[i][j] = w[0] * rho[i][j];
                feq1[i][j] = w[1] * rho[i][j];
                feq2[i][j] = w[2] * rho[i][j];
                feq3[i][j] = w[3] * rho[i][j];
                feq4[i][j] = w[4] * rho[i][j];
                feq5[i][j] = w[5] * rho[i][j];
                feq6[i][j] = w[6] * rho[i][j];
                feq7[i][j] = w[7] * rho[i][j];
                feq8[i][j] = w[8] * rho[i][j];


                f0[i][j] = f0[i][j] * (1.0 - omega) + omega * feq0[i][j];
                f1[i][j] = f1[i][j] * (1.0 - omega) + omega * feq1[i][j];
                f2[i][j] = f2[i][j] * (1.0 - omega) + omega * feq2[i][j];
                f3[i][j] = f3[i][j] * (1.0 - omega) + omega * feq3[i][j];
                f4[i][j] = f4[i][j] * (1.0 - omega) + omega * feq4[i][j];
                f5[i][j] = f5[i][j] * (1.0 - omega) + omega * feq5[i][j];
                f6[i][j] = f6[i][j] * (1.0 - omega) + omega * feq6[i][j];
                f7[i][j] = f7[i][j] * (1.0 - omega) + omega * feq7[i][j];
                f8[i][j] = f8[i][j] * (1.0 - omega) + omega * feq8[i][j];
            }
        }

        // ===================          2. streaming            ====================
        f0_old = f0;
        f1_old = f1;
        f2_old = f2;
        f3_old = f3;
        f4_old = f4;
        f5_old = f5;
        f6_old = f6;
        f7_old = f7;
        f8_old = f8;
        //
        // for (int j = M - 1; j > 0; --j) {
        //     for (int i = 0; i < M - 1; ++i) {
        //         f2[i][j] = f2[i][j - 1];
        //         f6[i][j] = f6[i - 1][j + 1];
        //     }
        // }
        //
        // for (int j = M - 1; j > 0; --j) {
        //     for (int i = M - 1; i > 0; --i) {
        //         f1[i][j] = f1[i - 1][j];
        //         f5[i][j] = f5[i - 1][j - 1];
        //     }
        // }
        //
        //
        // for (int j = 0; j < M - 1; ++j) {
        //     for (int i = M - 1; i > 0; --i) {
        //         f4[i][j] = f4[i][j + 1];
        //         f8[i][j] = f8[i - 1][j + 1];
        //     }
        // }
        //
        //
        // for (int j = 0; j < M - 1; ++j) {
        //     for (int i = 0; i < M - 1; ++i) {
        //         f3[i][j] = f3[i + 1][j];
        //         f7[i][j] = f7[i + 1][j + 1];
        //     }
        // }
        //
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                f0[i][j] = f0_old[i][j];

                f1[i][j] = (i > 0) ? f1_old[i - 1][j] : f1_old[i][j];
                f2[i][j] = (j > 0) ? f2_old[i][j - 1] : f2_old[i][j];
                f3[i][j] = (i < M - 1) ? f3_old[i + 1][j] : f3_old[i][j];
                f4[i][j] = (j < M - 1) ? f4_old[i][j + 1] : f4_old[i][j];

                f5[i][j] = (i > 0 && j > 0) ? f5_old[i - 1][j - 1] : f5_old[i][j];
                f6[i][j] = (i < M - 1 && j > 0) ? f6_old[i + 1][j - 1] : f6_old[i][j];
                f7[i][j] = (i < M - 1 && j < M - 1) ? f7_old[i + 1][j + 1] : f7_old[i][j];
                f8[i][j] = (i > 0 && j < M - 1) ? f8_old[i - 1][j + 1] : f8_old[i][j];
            }
        }

        // ===================          3. Boundary Conditions  ====================

        for (int j = 0; j < M; ++j) {
            f1[0][j] = w[1] * tw + w[3] * tw - f3[0][j];
            f5[0][j] = w[5] * tw + w[7] * tw - f7[0][j];
            f8[0][j] = w[8] * tw + w[6] * tw - f6[0][j];

            f3[M - 1][j] = -f1[M - 1][j];
            f6[M - 1][j] = -f8[M - 1][j];
            f7[M - 1][j] = -f5[M - 1][j];
        }
        for (int i = 0; i < M; ++i) {
            f4[i][M - 1] = -f2[i][M - 1];
            f7[i][M - 1] = -f5[i][M - 1];
            f8[i][M - 1] = -f6[i][M - 1];

            f1[i][0] = f1[i][1];
            f2[i][0] = f2[i][1];
            f3[i][0] = f3[i][1];
            f4[i][0] = f4[i][1];
            f5[i][0] = f5[i][1];
            f6[i][0] = f6[i][1];
            f7[i][0] = f7[i][1];
            f8[i][0] = f8[i][1];
        }
    }

    // ===================          4. Calculating Macroscopic quantities   ====================


    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {
            rho[i][j] = f0[i][j] + f1[i][j] + f2[i][j] + f3[i][j] + f4[i][j] + f5[i][j] + f6[i][j] + f7[i][j] + f8[
                            i][j];
        }
    }

    std::fstream file("data_d2q9.dat", std::ios::out);
    file << 'x' << ' ' << 'y' << ' ' << "rho" << std::endl;


    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {
            file << i << ' ' << j << ' ' << rho[i][j] << std::endl;
            // std::cout << rho[i][j];
        }
        // std::cout << std::endl;
    }
    file.close();
}
