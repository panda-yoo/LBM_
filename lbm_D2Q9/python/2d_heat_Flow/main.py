import numpy as np
import matplotlib.pyplot as plt
import constants as cst
from tqdm import tqdm

nx = 100
ny = 100
OMEGA = .4

nitr = 100


#  indices -> 0 and nx-1/ny-1 will be treated as buffer nodes
def cal_feq(feq, rho):
    for jy in range(ny):
        for ix in range(nx):
            for k in range(cst.NO_Q):
                feq[jy, ix, k] = cst.w[k] * rho[jy, ix]


def sim_2d_prob():
    # defining and initializing distributions/rho
    feq = np.zeros((ny, nx, cst.NO_Q))
    f = np.zeros((ny, nx, cst.NO_Q))
    rho = np.zeros((ny, nx))
    rho[ny // 2, nx // 2] = 1.0
    cal_feq(feq, rho)
    # plt.imshow(rho)
    # plt.show()
    for t in tqdm(range(nitr)):

        # collision
        for jy in range(1, ny - 1):
            for ix in range(1, nx - 1):
                f[jy, ix] = f[jy, ix] * (1.0 - OMEGA) + feq[jy, ix] * OMEGA

        # streaming
        f_temp = np.zeros_like(f)
        # f_temp[jy, ix, 1]
        # getting buffer nodes ready for streaming
        # left
        f[:, 0, 1] = f[:, nx - 2, 1]
        f[:, 0, 5] = f[:, nx - 2, 5]
        f[:, 0, 8] = f[:, nx - 2, 8]

        # right
        f[:, nx - 1, 3] = f[:, 1, 3]
        f[:, nx - 1, 6] = f[:, 1, 6]
        f[:, nx - 1, 7] = f[:, 1, 7]
        # down
        f[0, :, 8] = f[ny - 2, :, 8]
        f[0, :, 4] = f[ny - 2, :, 4]
        f[0, :, 7] = f[ny - 2, :, 7]
        # up
        f[ny - 1, :, 6] = f[1, :, 6]
        f[ny - 1, :, 2] = f[1, :, 2]
        f[ny - 1, :, 5] = f[1, :, 5]

        # corners
        f[0, 0, 5] = f[ny - 2, nx - 2, 5]
        f[0, nx - 1, 6] = f[ny - 2, 1, 6]
        f[ny - 1, 0, 8] = f[1, nx - 2, 8]
        f[ny - 1, nx - 1, 7] = f[1, 1, 7]

        for jy in range(1, ny - 1):
            for ix in range(1, nx - 1):
                f_temp[jy, ix, 1] = f[jy, ix - 1, 1]  # f1
                f_temp[jy, ix, 2] = f[jy - 1, ix, 2]  # f2
                f_temp[jy, ix, 3] = f[jy, ix + 1, 3]  # f3
                f_temp[jy, ix, 4] = f[jy + 1, ix, 4]  # f4

                f_temp[jy, ix, 5] = f[jy - 1, ix - 1, 5]  # f5
                f_temp[jy, ix, 6] = f[jy - 1, ix + 1, 6]  # f6
                f_temp[jy, ix, 7] = f[jy + 1, ix + 1, 7]  # f7
                f_temp[jy, ix, 8] = f[jy + 1, ix - 1, 8]  # f8

        f[:] = f_temp
        for jy in range(ny):
            for ix in range(nx):
                rho[jy, ix] = np.sum(f[jy, ix, :])

        cal_feq(feq, rho)

    fig, ax = plt.subplots()

    # ignoring buffer nodes while plotting
    ax.plot(np.arange(0, nx - 2, 1), rho[ny // 2, 1:-1])
    fig.savefig('results/1dplot.png')

    ax.imshow(rho[1:-1, 1:-1])
    fig.savefig('results/2dplot.png')

    plt.imshow(rho)
    plt.show()


if __name__ == "__main__":
    sim_2d_prob()
    pass
