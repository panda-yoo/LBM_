//
// Created by PAndA on 25-07-2025.
//
#include <fstream>
#include <array>
#include <iostream>
#include "../include/D1Q3.hpp"


void d1q3_sim() {
    std::array<double, N> f0;
    std::array<double, N> f1;
    std::array<double, N> f2;

    std::array<double, N> feq0;
    std::array<double, N> feq1;
    std::array<double, N> feq2;

    std::array<double, N> rho;

    f0.fill(0.0);
    f1.fill(0.0);
    f2.fill(0.0);

    float w0 = 4.0 / 6.0;
    float w1 = 1.0 / 6.0;
    float w2 = 1.0 / 6.0;

    float csq = 1.0 / 3.0;

    float alpha = 0.25;
    double omega = 1.0 / (3. * alpha + 0.5);

    // const float omega = 0.3;
    int no_times = 100;

    rho.fill(0.0);

    // adding disturbance
    // rho[N / 2] = 2.0;

    std::fstream file("data.dat", std::ios::out);

    file << " t " << " X " << " rho" << std::endl;

    // initializing the domain at i.e T(X= 0) = 1
    double rho_in = 1.0;
    f0[0] = w0 * rho_in;
    f1[0] = w1 * rho_in;
    f2[0] = w2 * rho_in;

    for (int time = 0; time < no_times; ++time) {
        // 1 ---- Collision

        for (int i = 0; i < N; ++i) {
            feq0[i] = w0 * rho[i];
            feq1[i] = w1 * rho[i];
            feq2[i] = w2 * rho[i];

            f0[i] = (1.0 - omega) * f0[i] + feq0[i] * omega;
            f1[i] = (1.0 - omega) * f1[i] + feq1[i] * omega;
            f2[i] = (1.0 - omega) * f2[i] + feq2[i] * omega;
        }


        // 2 ---- Streaming
        for (int i = 1; i <= N - 1; ++i) {
            f1[N - i] = f1[N - i - 1];
            f2[i - 1] = f2[i];
        }


        // Boundary conditions

        f1[0] = f2[0];
        f1[N - 1] = f1[N - 2];
        f2[N - 1] = f2[N - 2];


        // calculating macroscopic quantities
        for (int i = 0; i < N; ++i) {
            rho[i] = f0[i] + f1[i] + f2[i];
        }


        for (int i = 0; i < N; ++i) {
            if (file.is_open()) {
                file << time << ' ' << static_cast<double>(i) / (N - 1) << ' ' << rho[i] << std::endl;
            }
        }
    }


    file.close();
    std::cout << "file is closed " << std::endl;
}
