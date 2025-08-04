//
// Created by PAndA on 25-07-2025.
//

#ifndef D2Q9_HPP
#define D2Q9_HPP
#include <array>
constexpr int M = 50;

const std::array<double, 9> w{
    4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

void d2q9_sim();


#endif //D2Q9_HPP
