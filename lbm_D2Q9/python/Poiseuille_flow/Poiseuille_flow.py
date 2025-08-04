import numpy as np
import constants as cst
from typing import TypeAlias, Annotated
from numpy.typing import NDArray
import numba

NX = 50
NY = 50
dpdx = 1e-3
delta = 0.5
H = NY - 1.0 - 2.0 * delta
# buffer_ny = NY + 2
# distributionQ9: TypeAlias = Annotated[NDArray[np.float64], (cst.NO_Q, buffer_ny, NX)]
# plana2d: TypeAlias = Annotated[NDArray[np.float64], (buffer_ny, NX)]


distributionQ9: TypeAlias = Annotated[NDArray[np.float64], (cst.NO_Q, NY, NX)]
plana2d: TypeAlias = Annotated[NDArray[np.float64], (NY, NX)]


def is_wall(jy: int, ix: int):
    if jy in [0, NY]:
        return 1
    return 0


@numba.njit()
def collision(f: distributionQ9, feq: distributionQ9, ux: plana2d, uy: plana2d):
    for j in range(0, NY):
        for i in range(0, NX):
            for k in range(cst.NO_Q):
                c_dot_u = cst.cx[k] * ux[j][i] + cst.cy[k] * uy[j][i]
                term_gou_X = (1 - 1 / (2 * cst.TAU)) * cst.w[k] * (
                        (cst.cx[k] - ux[j][i]) / cst.CS_2 + (c_dot_u * cst.cx[k]) / cst.CS_2 ** 2) * dpdx
                f[k][j][i] = f[k][j][i] * (1.0 - cst.OMEGA) + cst.OMEGA * feq[k][j][i] + delta * term_gou_X


@numba.njit()
def streaming(f: distributionQ9):
    f_temp = f.copy()

    for j in range(0, NY):
        for i in range(0, NX):

            if i == NX - 1:
                # east‐pointing
                f_temp[3][j][i] = f[3][j][0]
                f_temp[6][j][i] = f[6][j - 1][0] if j > 0 else f[6][NY - 1][0]
                f_temp[7][j][i] = f[7][j + 1][0] if j < NY - 1 else f[7][0][0]

                # ——— WEST PERIODIC ———
            elif i == 0:
                # west‐pointing
                f_temp[1][j][i] = f[1][j][NX - 1]
                f_temp[5][j][i] = f[5][j - 1][NX - 1] if j > 0 else f[5][NY - 1][NX - 1]
                f_temp[8][j][i] = f[8][j + 1][NX - 1] if j < NY - 1 else f[8][0][NX - 1]

                # ——— SOUTH HALF-WAY BOUNCE-BACK ———
            elif j == 0:
                f_temp[2][j][i] = f[4][j][i]
                f_temp[5][j][i] = f[7][j][i]
                f_temp[6][j][i] = f[8][j][i]

                # ——— NORTH HALF-WAY BOUNCE-BACK ———
            elif j == NY - 1:
                f_temp[4][j][i] = f[2][j][i]
                f_temp[7][j][i] = f[5][j][i]
                f_temp[8][j][i] = f[6][j][i]
            if i == 0 and j == 0:
                pass
            # # ===============================================     EAST(periodic)       ================================================================
            # if i == NX - 1:
            #     f_temp[6][j][i] = f[6][j][0]
            #     f_temp[3][j][i] = f[3][j][0]
            #     f_temp[7][j][i] = f[7][j][0]
            #
            # # ===============================================     WEST(periodic)       ================================================================
            # elif i == 0:
            #     f_temp[5][j][i] = f[5][j][NX - 1]
            #     f_temp[1][j][i] = f[1][j][NX - 1]
            #     f_temp[8][j][i] = f[8][j][NX - 1]
            #
            # # ===============================================     SOUTH(halfway bounceback)       ================================================================
            # elif j == 0:
            #     f_temp[2][j][i] = f[4][j][i]
            #     f_temp[5][j][i] = f[7][j][i]
            #     f_temp[6][j][i] = f[8][j][i]
            #
            # # ===============================================     NORTH(halfway bounceback)       ================================================================
            # elif j == NY - 1:
            #     f_temp[7][j][i] = f[5][j][i]
            #     f_temp[4][j][i] = f[2][j][i]
            #     f_temp[8][j][i] = f[6][j][i]

            else:
                f_temp[1][j][i] = f[1][j][i - 1]
                f_temp[2][j][i] = f[2][j - 1][i]
                f_temp[3][j][i] = f[3][j][i + 1]
                f_temp[4][j][i] = f[4][j + 1][i]

                f_temp[5][j][i] = f[5][j - 1][i - 1]
                f_temp[6][j][i] = f[6][j - 1][i + 1]
                f_temp[7][j][i] = f[7][j + 1][i + 1]
                f_temp[8][j][i] = f[8][j + 1][i - 1]

    return f_temp


@numba.njit()
def calculate_f_eq(feq: distributionQ9, rho: plana2d, ux: plana2d, uy: plana2d):
    for j in range(0, NY):
        for i in range(0, NX):

            u_dot_u = ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i]
            for k in range(cst.NO_Q):
                u_dot_c = ux[j][i] * cst.cx[k] + uy[j][i] * cst.cy[k]

                feq[k][j][i] = (
                        cst.w[k]
                        * rho[j][i]
                        * (
                                1.0
                                + u_dot_c / cst.CS_2
                                + (u_dot_c * u_dot_c) / (2 * cst.CS_2 * cst.CS_2)
                                - u_dot_u / (2 * cst.CS_2)
                        )
                )


@numba.njit()
def calculate_macroscopic_quan(f: distributionQ9, rho: plana2d, ux: plana2d, uy: plana2d):
    for j in range(0, NY):
        for i in range(0, NX):

            sum_rho = 0.0
            sum_ux = 0.0
            sum_uy = 0.0
            for k in range(cst.NO_Q):
                sum_rho = sum_rho + f[k][j][i]
                sum_ux = sum_ux + f[k][j][i] * cst.cx[k] + dpdx / 2.0
                sum_uy = sum_uy + f[k][j][i] * cst.cy[k]

            if sum_rho < 1e-5:
                sum_rho = 1e-2
            rho[j][i] = sum_rho

            ux[j][i] = sum_ux / sum_rho
            uy[j][i] = sum_uy / sum_rho

            # if ux[j][i] > 1.0:
            #     ux[j][i] = 1.0
            # if ux[j][i] < -1.0:
            #     ux[j][i] = -1.0
            #
            # if uy[j][i] > 1.0:
            #     uy[j][i] = 1.0
            # if uy[j][i] < -1.0:
            #     uy[j][i] = -1.0
            # ux[j][i] = np.clip(sum_ux / sum_rho, -1.0, 1.0)
            # uy[j][i] = np.clip(sum_uy / sum_rho, -1.0, 1.0)

            # ux[j][i] = sum_ux / sum_rho
            # uy[j][i] = sum_uy / sum_rho
