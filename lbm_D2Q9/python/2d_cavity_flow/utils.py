from typing import Annotated, TypeAlias
import numpy as np
from numpy.typing import NDArray
import parameters as para

plane2d: TypeAlias = Annotated[NDArray[np.float64], (para.ny1, para.nx1)]

distributionQ9: TypeAlias = Annotated[NDArray[np.float64], (para.NO_Q, para.ny1, para.nx1)]


def Init_Eq(f: distributionQ9, rho: plane2d, ux: plane2d, uy: plane2d):
    for j in range(para.ny1):
        for i in range(para.nx1):
            rho[j][i] = para.rho0
            ux[j][i] = para.ux0
            uy[j][i] = para.uy0
            for k in range(para.NO_Q):
                f[k][j][i] = feq(rho=rho[j][i], ux=ux[j][i], uy=uy[j][i], k=k)

def feq(rho: np.float64, ux: np.float64, uy: np.float64, k: int):
    c_dot_u = para.cx[k] * ux + para.cy[k] * uy
    u_dot_u = ux * ux + uy * uy

    return para.w[k] * rho * (1 + 3 * c_dot_u + 4.5 * c_dot_u * c_dot_u - 1.5 * u_dot_u)


def coll_BKG(f: distributionQ9, rho: np.float64, ux: np.float64, uy: np.float64, k: int):
    f_post = np.zeros((para.NO_Q, para.nx1, para.ny1))
    for j in range(para.ny1):
        for i in range(para.nx1):

            for k in range(para.NO_Q):
                FEQ = feq(rho=rho[j][i], ux=ux[j][i], uy=uy[j][i], k=k)

                f_post[j][i][k] = f[j][i][k] - (f[j][i][k] - FEQ) / para.TAU


def streaming(f: distributionQ9, f_post: distributionQ9):
    for j in range(para.ny1):
        for i in range(para.nx1):
            for k in range(para.NO_Q):
                jd_ = j - para.cy[k]
                id_ = i - para.cx[k]
                if 0 <= jd_ <= para.NY and 0 <= id_ <= para.NX:
                    f[j][i][k] = f_post[jd_][id_][k]


def bounce_back(f: distributionQ9, f_post: distributionQ9, rho: plane2d, boundary: plane2d):
    # j = Ny: top plate
    for i in range(para.nx1):
        f[para.NY][i][4] = f_post[para.NY][i][2]
        f[para.NY][i][7] = f_post[para.NY][i][5] + 6 * rho[para.NY][i] * para.w[7] * para.cx[7] * para.uw
        f[para.NY][i][8] = f_post[para.NY][i][6] + 6 * rho[para.NY][i] * para.w[8] * para.cx[8] * para.uw

    #  j = 0: bottom plate
    for i in range(para.nx1):
        f[0][i][2] = f_post[0][i][4]
        f[0][i][5] = f_post[0][i][7]
        f[0][i][6] = f_post[0][i][8]

    #  i = 0: left wall
    for j in range(para.ny1):
        f[j][0][1] = f_post[j][0][3]
        f[j][0][5] = f_post[j][0][7]
        f[j][0][8] = f_post[j][0][6]

    #  i = Nx: right wall
    for j in range(para.ny1):
        f[j][para.NX][3] = f_post[j][para.NX][1]
        f[j][para.NX][7] = f_post[j][para.NX][5]
        f[j][para.NX][6] = f_post[j][para.NX][8]


def den_Vel(f: distributionQ9, rho: plane2d, ux: plane2d, uy: plane2d):
    for j in range(para.ny1):
        for i in range(para.nx1):
            rho[j][i] = f[j][i][0] + f[j][i][1] + f[j][i][2] + f[j][i][3] + f[j][i][4] + f[j][i][5] + f[j][i][6] + \
                        f[j][i][7] + f[j][i][8]
            ux[j][i] = (f[j][i][1] + f[j][i][5] + f[j][i][8] - f[j][i][3] - f[j][i][6] - f[j][i][7]) / rho[j][i]

            uy[j][i] = (f[j][i][5] + f[j][i][6] + f[j][i][2] - f[j][i][7] - f[j][i][8] - f[j][i][4]) / rho[j][i]
