//
// Created by PAndA on 26-07-2025.
//
#include <array>

#ifndef D2Q9_CYLINDER_HPP
#define D2Q9_CYLINDER_HPP


const unsigned int nx = 50;
const unsigned int ny = 20;
const unsigned int q = 9;


typedef std::array<std::array<double, nx>, ny> plane2d;
typedef std::array<std::array<std::array<double, nx>, ny>, q> distributionsQ9;


const double centre[2] = {(nx - 1.0) / 4, (ny - 1.0) / 2};

const double cs2 = 1.0 / 3.0;
const double rad = 4.0;
const double VELOCITY_X = 0.01;
const double OMEGA = 1.0;
const std::array<double, q> weight{
    4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};
const int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

const std::array<int, q> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
const std::array<int, q> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};


double distance(int i, int j, double n1, double n2);

void cylinder_create(plane2d &cyl);


void show(plane2d &cyl);

void Calculating_u_Comp(const plane2d &cylinder,
                        plane2d &rho,
                        distributionsQ9 &f,
                        plane2d &ux,
                        plane2d &uy);

void CalculatingMacroscopicQuantities(const plane2d &cylinder,
                                      plane2d &rho,
                                      distributionsQ9 &f);

void collision(const plane2d &cylinder,
               distributionsQ9 &f,
               distributionsQ9 &feq);


void streaming(const plane2d &cylinder,
               distributionsQ9 &f);

void f_equilibrium(const plane2d &cylinder,
                   const plane2d &rho,
                   distributionsQ9 &feq,
                   plane2d &ux,
                   plane2d &uy);


void boundary(const plane2d &cylinder,
              distributionsQ9 &f,
              const plane2d &rho);

void sim_Cylinder_d2q9();


#endif //D2Q9_CYLINDER_HPP
