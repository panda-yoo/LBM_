import numpy as np
from typing import Annotated, TypeAlias
from numpy.typing import NDArray
import constants as cst

# defining types for simplicity/readability (with including buffer nodes)

distribution_Q9: TypeAlias = Annotated[
    NDArray[np.float64], (cst.NO_Q, cst.buffer_ny, cst.buffer_nx)
]
plane2d: TypeAlias = Annotated[NDArray[np.float64], (cst.buffer_ny, cst.buffer_nx)]

plane2d_cylinder: TypeAlias = Annotated[
    NDArray[np.int8], (cst.buffer_ny, cst.buffer_nx)
]


def distance(i: int, j: int, nx: int, ny: int):
    return ((i - nx) ** 2 + (j - ny) ** 2) ** 0.5


def generate_isCylinder(iscylinde: plane2d_cylinder):
    for j in range(cst.buffer_ny):
        for i in range(cst.buffer_nx):
            if (
                    distance(i=i, j=j, nx=cst.CENTER_CYLINDER[0], ny=cst.CENTER_CYLINDER[1])
                    < cst.CYLINDER_RADIUS
            ):
                iscylinde[j][i] = 1


def streaming(iscylinde: plane2d_cylinder, f: distribution_Q9) -> distribution_Q9:
    # excluding buffer nodes
    f_temp = f.copy()
    for j in range(1, cst.NY + 1):
        for i in range(1, cst.NX + 1):
            if iscylinde[j][i] == 0:
                f_temp[1][j][i] = f[1][j][i - 1]
                f_temp[2][j][i] = f[2][j - 1][i]
                f_temp[3][j][i] = f[3][j][i + 1]
                f_temp[4][j][i] = f[4][j + 1][i]

                f_temp[5][j][i] = f[5][j - 1][i - 1]
                f_temp[6][j][i] = f[6][j - 1][i + 1]
                f_temp[7][j][i] = f[7][j + 1][i + 1]
                f_temp[8][j][i] = f[8][j + 1][i - 1]

    return f_temp


def collision(iscylinde: plane2d_cylinder, f: distribution_Q9, feq: distribution_Q9):
    for j in range(1, cst.NY + 1):
        for i in range(1, cst.NX + 1):
            if iscylinde[j][i] == 0:
                for k in range(cst.NO_Q):
                    f[k][j][i] = (1.0 - cst.OMEGA) * f[k][j][i] + cst.OMEGA * feq[k][j][
                        i
                    ]


def calculate_f_equilibrium(
        iscylinde: plane2d_cylinder,
        feq: distribution_Q9,
        rho: plane2d,
        ux: plane2d,
        uy: plane2d,
):
    for j in range(1, cst.NY + 1):
        for i in range(1, cst.NX + 1):
            if iscylinde[j][i] == 0:
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


def calculate_microsopic_quantities(
        iscylinde: plane2d_cylinder,
        f: distribution_Q9,
        rho: plane2d,
        ux: plane2d,
        uy: plane2d,
):
    for j in range(1, cst.NY + 1):
        for i in range(1, cst.NX + 1):
            if iscylinde[j][i] == 0:
                sum_rho = 0.0
                sum_ux = 0.0
                sum_uy = 0.0
                for k in range(cst.NO_Q):
                    sum_rho = sum_rho + f[k][j][i]
                    sum_ux = sum_ux + f[k][j][i] * cst.cx[k]
                    sum_uy = sum_uy + f[k][j][i] * cst.cy[k]
                # if sum_rho < 1e-5:
                #     sum_rho = 1e-2
                if not np.isfinite(sum_rho) or sum_rho < 1e-5:
                    sum_rho = 1e-1
                rho[j][i] = sum_rho

                ux[j][i] = np.clip(sum_ux / sum_rho, -1.0, 1.0)
                uy[j][i] = np.clip(sum_uy / sum_rho, -1.0, 1.0)

                # ux[j][i] = sum_ux / sum_rho
                # uy[j][i] = sum_uy / sum_rho

    # preparing ghost nodes for streaming


def apply_zou_he_boundary(f: distribution_Q9, rho: plane2d):
    for j in range(0, cst.buffer_ny):
        # ===========================================       WEST BOUNDARY       =====================================================
        rho_temp = (
                           f[0][j][1]
                           + f[2][j][1]
                           + f[4][j][1]
                           + 2 * (f[3][j][1] + f[6][j][1] + f[7][j][1])
                   ) / (1 - cst.INFLOW_VELOCITY_UX)

        f[1][j][0] = f[3][j][1] + (2.0 / 3.0) * rho_temp * cst.INFLOW_VELOCITY_UX

        f[5][j][0] = (
                f[7][j][1]
                - (1.0 / 2.0) * (f[2][j][1] - f[4][j][1])
                + (1.0 / 6.0) * rho_temp * cst.INFLOW_VELOCITY_UX
        )

        f[8][j][0] = (
                f[6][j][1]
                + (1.0 / 2.0) * (f[2][j][1] - f[4][j][1])
                + (1.0 / 6.0) * rho_temp * cst.INFLOW_VELOCITY_UX
        )
        # ===========================================       EAST BOUNDARY       =====================================================
        f[6][j][cst.buffer_nx - 1] = f[6][j][cst.buffer_nx - 2]
        f[3][j][cst.buffer_nx - 1] = f[3][j][cst.buffer_nx - 2]
        f[7][j][cst.buffer_nx - 1] = f[7][j][cst.buffer_nx - 2]


def apply_wall_bounceback(f: distribution_Q9, rho: plane2d):
    for i in range(0, cst.buffer_nx):
        # ===========================================       SOUTH BOUNDARY       =====================================================
        f[2][0][i] = f[4][1][i]
        f[6][0][i] = f[8][1][i]
        f[5][0][i] = f[7][1][i]

        # ===========================================       NORTH BOUNDARY       =====================================================
        f[7][cst.buffer_ny - 1][i] = f[5][cst.buffer_ny - 2][i]
        f[4][cst.buffer_ny - 1][i] = f[2][cst.buffer_ny - 2][i]
        f[8][cst.buffer_ny - 1][i] = f[6][cst.buffer_ny - 2][i]


# def apply_obstacle_bounceback(
#     iscylinder: plane2d_cylinder, f: distribution_Q9, rho: plane2d
# ):
#     for j in range(1, cst.NY + 1):
#         for i in range(1, cst.NX + 1):
#             if iscylinder[j][i] == 1:
#                 f[3][j][i] = f[1][j][i - 1]
#                 f[4][j][i] = f[2][j - 1][i]
#                 f[1][j][i] = f[3][j][i + 1]
#                 f[2][j][i] = f[4][j + 1][i]


#                 f[7][j][i] = f[5][j - 1][i - 1]
#                 f[8][j][i] = f[6][j - 1][i + 1]
#                 f[5][j][i] = f[7][j + 1][i + 1]
#                 f[6][j][i] = f[8][j + 1][i - 1]
def apply_obstacle_bounceback(
        iscylinder: plane2d_cylinder, f: distribution_Q9, rho: plane2d
):
    opposite = [0, 3, 4, 1, 2, 7, 8, 5, 6]
    for j in range(1, cst.NY + 1):
        for i in range(1, cst.NX + 1):
            if iscylinder[j][i] == 1:
                for k in range(1, cst.NO_Q):
                    f[k][j][i] = f[opposite[k]][j][i]

# def boundary():
#     pass
