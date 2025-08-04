//
// Created by PAndA on 26-07-2025.
//

#include "../include/D2Q9_cylinder.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>


// typedef std::array<std::array<double, ny>, nx> plane2d;
// typedef std::array<std::array<double, nx>, ny> plane2d;
// typedef std::array<std::array<std::array<double, nx>, ny>, q> distributionsQ9;

double distance(int i, int j, double n1, double n2) {
    double d2 = std::pow((n1 - i), 2) + std::pow(n2 - j, 2);
    return std::pow(d2, 0.5);
}

void cylinder_create(plane2d &cyl) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (distance(i, j, centre[0], centre[1]) < rad) {
                cyl[j][i] = 1.0;
            } else {
                cyl[j][i] = 0.0;
            }
        }
    }
}

void show(plane2d &cyl) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            std::cout << cyl[j][i] << ' ';
        }
        std::cout << std::endl;
    }
}

void Calculating_u_Comp(const plane2d &cylinder,
                        plane2d &rho,
                        distributionsQ9 &f,
                        plane2d &ux,
                        plane2d &uy) {
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            // 1 for cylinder is there
            // 0 for cylinder is not there
            if (cylinder[j][i] == 0.0) {
                double sumx = 0.0;
                double sumy = 0.0;

                for (int k = 0; k < q; ++k) {
                    sumx = sumx + f[k][j][i] * cx[k];
                    sumy = sumy + f[k][j][i] * cy[k];
                }
                if (rho[j][i] < 1e-7) {
                    rho[j][i] = 1e-5;
                }
                ux[j][i] = sumx / rho[j][i];
                uy[j][i] = sumy / rho[j][i];
            } else {
                ux[j][i] = 0.0;
                uy[j][i] = 0.0;
            }
        }
    }
}

void CalculatingMacroscopicQuantities(const plane2d &cylinder,
                                      plane2d &rho,
                                      distributionsQ9 &f) {
    double sum = 0.0;
    for (int j = 1; j < ny - 1; ++j) {
        //WEST side
        // rho[0][j] = (f[0][0][j] + f[2][0][j] + f[4][0][j]
        //              + 2 * (f[3][0][j] + f[6][0][j] + f[7][0][j])) / (1.0 - ux);

        rho[j][0] = (f[0][j][0] + f[2][j][0] + f[4][j][0]
                     + 2 * (f[3][j][0] + f[6][j][0] + f[7][j][0])) / (1.0 - VELOCITY_X);
    }

    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            // 1 for cylinder is there
            // 0 for cylinder is not there
            if (cylinder[j][i] == 0.0) {
                sum = 0.0;
                for (int k = 0; k < q; ++k) {
                    sum = sum + f[k][j][i];
                }
                rho[j][i] = sum;
            } else {
                rho[j][i] = 0.0;
            }
        }
    }

    // =============================================================================================================================
    for (int j = 0; j < ny; ++j) {
        // Iterate over all rows at the inlet
        // Check if it's not a corner where other BCs might apply first
        if (cylinder[j][0] == 0.0) {
            // Only if it's a fluid node at the inlet
            // This is a common form for inlet density reconstruction where ux is fixed to VELOCITY_X
            rho[j][0] = (f[0][j][0] + f[2][j][0] + f[4][j][0] + 2.0 * (f[3][j][0] + f[6][j][0] + f[7][j][0])) / (
                            1.0 - VELOCITY_X);
            // NOTE: This assumes uy = 0 and fixed ux.
            // Ensure this formula is correct for your chosen inlet BC method.
        }
        // =============================================================================================================================
    }
}

void collision(const plane2d &cylinder,
               distributionsQ9 &f,
               distributionsQ9 &feq) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (cylinder[j][i] == 0.0) {
                for (int k = 0; k < q; ++k) {
                    f[k][j][i] = f[k][j][i] * (1.0 - OMEGA) + OMEGA * feq[k][j][i];
                }
            }
        }
    }
}

//
// void streaming(const plane2d &cylinder,
//                distributionsQ9 &f) {
//     distributionsQ9 f_temp;
//
//     for (int j = 1; j < ny - 1; ++j) {
//         for (int i = 1; i < nx - 1; ++i) {
//             if (cylinder[j][i] == 0.0) {
//                 f_temp[0][j][i] = f[0][j][i];
//
//                 f_temp[1][j][i] = f[1][j][i - 1];
//
//                 f_temp[2][j][i] = f[2][j - 1][i];
//
//                 f_temp[3][j][i] = f[3][j][i + 1];
//
//                 f_temp[4][j][i] = f[4][j + 1][i];
//
//                 f_temp[5][j][i] = f[5][j - 1][i - 1];
//
//                 f_temp[6][j][i] = f[6][j - 1][i + 1];
//
//                 f_temp[7][j][i] = f[7][j + 1][i + 1];
//
//                 f_temp[8][j][i] = f[8][j + 1][i - 1];
//             } else {
//                 for (int k = 0; k < q; ++k) {
//                     f_temp[k][j][i] = 0.0;
//                 }
//             }
//         }
//     }
//     std::swap(f, f_temp);
// }
// D2Q9_cylinder.cpp snippet (streaming function, "hardcoded" for boundary cases)

void streaming(const plane2d &cylinder,
               distributionsQ9 &f) {
    // Create a temporary array for the *next* state after streaming.
    // This is essential to avoid in-place update issues during streaming.
    distributionsQ9 f_temp;

    // Loop over ALL target cells in the domain (i, j)
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // For each direction k
            for (int k = 0; k < q; ++k) {
                // Calculate the source coordinates for population f_k ending at (j, i)
                int prev_i = i - cx[k];
                int prev_j = j - cy[k];

                // Initialize the temporary distribution for this (k, j, i) to 0.0.
                // This means, by default, no population is streamed into this slot
                // unless it comes from a valid fluid source.
                // Populations coming from outside boundaries or solid nodes will
                // be explicitly set by the 'boundary' function later.
                f_temp[k][j][i] = 0.0;

                // Check if the source cell (prev_j, prev_i) is within the main grid boundaries
                if (prev_i >= 0 && prev_i < nx && prev_j >= 0 && prev_j < ny) {
                    // Check if the source cell is a fluid node (i.e., not part of the cylinder)
                    // Populations from solid nodes (cylinder) do not stream *out* into fluid cells.
                    // Instead, incoming populations to solid nodes are handled by bounce-back.
                    if (cylinder[prev_j][prev_i] == 0.0) {
                        // If source is a fluid node
                        f_temp[k][j][i] = f[k][prev_j][prev_i]; // Stream the population from the fluid source
                    }
                    // If the source (prev_j, prev_i) is a solid node (cylinder[prev_j][prev_i] == 1.0),
                    // then f_temp[k][j][i] remains 0.0 as initialized. The actual
                    // value for this population will be determined by the 'boundary' function later.
                }
                // If the source (prev_j, prev_i) is outside the domain boundaries,
                // then f_temp[k][j][i] remains 0.0 as initialized. The actual
                // value for this population will be determined by the 'boundary' function later.
            }
        }
    }
    // CRITICAL FIX: Swap f and f_temp *once* after all streaming for all cells is complete.
    // This moves the newly streamed distributions into 'f' for the next step.
    std::swap(f, f_temp);
}


void f_equilibrium(const plane2d &cylinder,
                   const plane2d &rho,
                   distributionsQ9 &feq,
                   plane2d &ux,
                   plane2d &uy) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (cylinder[j][i] == 0.0) {
                double u_dot_u = ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i];
                for (int k = 0; k < q; ++k) {
                    double u_dot_c = ux[j][i] * cx[k] + uy[j][i] * cy[k];

                    feq[k][j][i] = weight[k] * rho[j][i] * (1.0 + u_dot_c / cs2
                                                            + (u_dot_c * u_dot_c) / (2 * (cs2 * cs2))
                                                            - u_dot_u / (2 * cs2));
                }
            } else {
                // Solid node
                // For solid nodes, feq might not be directly relevant or set to 0.0,
                // depending on how bounce-back is implemented.
                for (int k = 0; k < q; ++k) {
                    feq[k][j][i] = 0.0; // Or keep existing values if not explicitly used
                }
            }
        }
    }
}

//
// void boundary(const plane2d &cylinder,
//               distributionsQ9 &f,
//               const plane2d &rho) {
//     for (int i = 0; i < nx ; ++i) {
//
//         //south side
//         f[6][0][i] = f[8][0][i];
//         f[2][0][i] = f[4][0][i];
//         f[5][0][i] = f[7][0][i];
//
//         //north side
//         f[7][ny - 1][i] = f[5][ny - 1][i];
//         f[4][ny - 1][i] = f[2][ny - 1][i];
//         f[8][ny - 1][i] = f[6][ny - 1][i];
//     }
//     for (int j = 1; j < ny - 1; ++j) {
//         //WEST side
//         f[1][j][0] = f[3][j][0] + (2.0 / 3.0) * rho[j][0] * VELOCITY_X;
//         f[5][j][0] = f[7][j][0] - (1.0 / 2.0) * (f[2][j][0] - f[4][j][0])
//                      + (1.0 / 6.0) * rho[j][0] * VELOCITY_X;
//         f[8][j][0] = f[6][j][0] + (1.0 / 2.0) * (f[2][j][0] - f[4][j][0])
//                      + (1.0 / 6.0) * rho[j][0] * VELOCITY_X;
//
//
//         //EAST side
//         f[6][j][nx - 1] = f[6][j][nx - 2];
//         f[7][j][nx - 1] = f[7][j][nx - 2];
//         f[3][j][nx - 1] = f[3][j][nx - 2];
//     }
//     //CYLINDER
//     for (int j = 1; j < ny - 1; ++j) {
//         for (int i = 1; i < nx - 1; ++i) {
//             if (cylinder[j][i] == 1.0) {
//                 std::array<double, q> temp;
//                 for (int k = 0; k < q; ++k)
//                     temp[k] = f[k][j][i];
//
//                 for (int k = 1; k < q; ++k)
//                     f[k][j][i] = temp[opposite[k]];
//             }
//         }
//     }
// }
// D2Q9_cylinder.cpp snippet (boundary)

void boundary(const plane2d &cylinder,
              distributionsQ9 &f,
              plane2d &rho) {
    // uy is needed for inlet BC


    // 1. South Wall (j=0, No-slip Bounce-back)
    for (int i = 0; i < nx; ++i) {
        if (cylinder[0][i] == 0.0) {
            // Only apply if it's a fluid node
            f[5][0][i] = f[7][0][i]; // f5 (SW) comes from f7 (NE)
            f[2][0][i] = f[4][0][i]; // f2 (N) comes from f4 (S)
            f[6][0][i] = f[8][0][i]; // f6 (NW) comes from f8 (SE)
        }
    }

    // 2. North Wall (j=ny-1, No-slip Bounce-back)
    for (int i = 0; i < nx; ++i) {
        if (cylinder[ny - 1][i] == 0.0) {
            // Only apply if it's a fluid node
            f[7][ny - 1][i] = f[5][ny - 1][i]; // f7 (SE) comes from f5 (NW)
            f[4][ny - 1][i] = f[2][ny - 1][i]; // f4 (S) comes from f2 (N)
            f[8][ny - 1][i] = f[6][ny - 1][i]; // f8 (NE) comes from f6 (SW)
        }
    }

    // 3. West Wall (i=0, Inlet - Guo et al. type with fixed velocity)
    for (int j = 0; j < ny; ++j) {
        if (cylinder[j][0] == 0.0) {
            // Only if it's a fluid node at inlet
            // First, calculate density at this inlet point (if not done in MacroscopicQuantities)
            // This is critical if the boundary method *defines* rho rather than uses a pre-calculated one.
            // If CalculatingMacroscopicQuantities already sets it, remove this.
            // This line fixes the typo and correctly reconstructs rho at the inlet based on fixed ux
            rho[j][0] = (f[0][j][0] + f[2][j][0] + f[4][j][0] + 2.0 * (f[3][j][0] + f[6][j][0] + f[7][j][0])) / (
                            1.0 - VELOCITY_X);

            // Reconstruct unknown populations (those streaming IN from west, i.e., directions 1, 5, 8)
            double ux_in = VELOCITY_X; // Fixed x-velocity at inlet
            double uy_in = 0.0; // Fixed y-velocity at inlet

            // Equilibrium parts (reused in multiple equations)
            double cu_eq_1 = cx[1] * ux_in + cy[1] * uy_in;
            double cu_eq_5 = cx[5] * ux_in + cy[5] * uy_in;
            double cu_eq_8 = cx[8] * ux_in + cy[8] * uy_in;
            double usq_eq = ux_in * ux_in + uy_in * uy_in;

            // These are based on Guo et al. for fixed velocity inlet (ux, uy, rho at boundary)
            f[1][j][0] = f[3][j][0] + (2.0 / 3.0) * rho[j][0] * ux_in;
            f[5][j][0] = f[7][j][0] - 0.5 * (f[2][j][0] - f[4][j][0]) + 0.5 * (
                             rho[j][0] * (cu_eq_5 / cs2 + 0.5 * cu_eq_5 * cu_eq_5 / (cs2 * cs2) - 0.5 * usq_eq / cs2));
            // Reconstruct f5
            f[8][j][0] = f[6][j][0] + 0.5 * (f[2][j][0] - f[4][j][0]) + 0.5 * (
                             rho[j][0] * (cu_eq_8 / cs2 + 0.5 * cu_eq_8 * cu_eq_8 / (cs2 * cs2) - 0.5 * usq_eq / cs2));
            // Reconstruct f8

            // Alternative simpler reconstruction for f5 and f8 (if uy=0)
            // f[5][j][0] = f[7][j][0] + 0.5 * (f[4][j][0] - f[2][j][0]) + (1.0/6.0) * rho[j][0] * ux_in;
            // f[8][j][0] = f[6][j][0] + 0.5 * (f[2][j][0] - f[4][j][0]) + (1.0/6.0) * rho[j][0] * ux_in;
            // The exact formulas vary slightly with specific implementations of Guo/Zou-He.
            // Double-check these against your chosen Guo method reference.
        }
    }

    // 4. East Wall (i=nx-1, Outlet - Zero-Gradient Extrapolation)
    for (int j = 0; j < ny; ++j) {
        if (cylinder[j][nx - 1] == 0.0) {
            // Only if it's a fluid node at outlet
            // Populations streaming IN from east are known (e.g., f[3], f[6], f[7])
            // Populations streaming OUT to east (f[1], f[5], f[8]) are unknown.
            // Extrapolate unknown populations from interior cells
            for (int k = 0; k < q; ++k) {
                // If this population (k) streams FROM the domain (cx[k] < 0 for East wall)
                // it is unknown, so extrapolate from left neighbor.
                if (cx[k] < 0) {
                    // Populations 3, 6, 7 stream from East
                    f[k][j][nx - 1] = f[k][j][nx - 2];
                }
            }
            // For populations streaming IN from west, (cx[k] > 0), those are known
            // and should already be there from streaming.
        }
    }

    // 5. Cylinder (Solid Obstacle - Bounce-back)
    // This part should be applied *after* streaming, to the distributions that have just streamed.
    // It should iterate over all cells, find solid ones, and apply bounce-back.
    // A common way is to apply bounce-back to the distributions for *fluid nodes adjacent to solid nodes*.
    // OR: for solid nodes directly, take incoming distributions from fluid neighbors and send back.
    // Your current approach is for solid nodes, which is fine if consistent.
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (cylinder[j][i] == 1.0) {
                // If it's a cylinder (solid) node
                // Store incoming distributions to a temporary array first,
                // as you might overwrite before reading all needed values.
                distributionsQ9 f_copy = f; // Make a local copy of f for this cell to get incoming values

                for (int k = 0; k < q; ++k) {
                    // This applies bounce-back to the populations at the solid node itself.
                    // This assumes that streaming has already happened,
                    // and 'f' now contains the populations after streaming.
                    f[k][j][i] = f_copy[opposite[k]][j][i];
                }
            }
        }
    }
}

void sim_Cylinder_d2q9() {
    const unsigned int ntimes = 1000;


    distributionsQ9 feq;

    plane2d cylinder;
    cylinder_create(cylinder);

    plane2d rho;


    std::array<double, nx> rows;

    rows.fill(1.0);

    rho.fill(rows);

    rows.fill(0.0);

    plane2d ux;
    plane2d uy;

    ux.fill(rows);
    uy.fill(rows);

    for (int j = 0; j < ny; ++j) {
        ux[j][0] = VELOCITY_X;
    }


    f_equilibrium(cylinder, rho, feq, ux, uy);

    distributionsQ9 f = feq;


    std::fstream file("data_vortex.dat", std::ios::out);

    // file << "t" << " " << "x" << " " << "y" << " " << "rho" << " " << "ux" << " " << "uy" << std::endl;
    file << "x" << " " << "y" << " " << "rho" << " " << "ux" << " " << "uy" << std::endl;

    for (int t = 0; t < ntimes; ++t) {
        f_equilibrium(cylinder, rho, feq, ux, uy);
        // =================================        1. Collision        =======================================
        collision(cylinder, f, feq);

        // =================================        2. Streaming        =======================================
        streaming(cylinder, f);

        // =================================        3. Boundary        =======================================

        boundary(cylinder, f, rho);
        // =================================        4. Calculating Macroscopic quantities        =======================================
        CalculatingMacroscopicQuantities(cylinder, rho, f);
        Calculating_u_Comp(cylinder, rho, f, ux, uy);
    }
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            // file << "t" << " " << "x" << " " << "y" << " " << "rho" << " " << "ux" << " " << "uy" << std::endl;

            // file << t << " " << i << " " << j << " " << rho[i][j] << " " << ux[i][j] << " " << uy[i][j] <<
            //         std::endl;
            file << i << " " << j << " " << rho[i][j] << " " << ux[i][j] << " " << uy[i][j] <<
                    std::endl;
        }
    }
    file.close();
}
