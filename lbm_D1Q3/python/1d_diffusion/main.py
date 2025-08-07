import matplotlib.pyplot as plt
import numpy as np


def sim():
    m = 101

    OMEGA = 0.5
    niter = 1000

    w0 = 4. / 6.
    w1 = 1. / 6.
    w2 = 1. / 6.
    cs2 = 1. / 3.
    # initialise f,rho
    f0 = np.zeros(m, dtype=np.float64)
    f1 = np.zeros(m, dtype=np.float64)
    f2 = np.zeros(m, dtype=np.float64)

    rho = np.zeros(m, dtype=np.float64)

    # feq_i = w_k * rho_i
    # initializing distributions
    for i in range(m):
        f0[i] = w0 * rho[i]
        f1[i] = w1 * rho[i]
        f2[i] = w2 * rho[i]
    for t in range(niter + 1):
        # collisions
        for i in range(m):
            f0eq = w0 * rho[i]
            f1eq = w1 * rho[i]

            f0[i] = f0[i] * (1.0 - OMEGA) + OMEGA * f0eq
            f1[i] = f1[i] * (1.0 - OMEGA) + OMEGA * f1eq
            f2[i] = f2[i] * (1.0 - OMEGA) + OMEGA * f1eq
        # streaming
        f1_temp = f1.copy()
        f2_temp = f2.copy()
        for i in range(1, m - 1):
            f1_temp[i] = f1[i - 1]
            f2_temp[i] = f2[i + 1]

        f1 = f1_temp
        f2 = f2_temp
        # boundary condition

        # its like maintaining constant value at boundary
        f1[0] = 1.0 - (f0[0] + f2[0])
        # just copying distribution f2[last_node-2] to f2[last_node-1]
        f2[m - 1] = f2[m - 2]
        for i in range(m):
            rho[i] = f0[i] + f1[i] + f2[i]

        if t % 100 == 0:
            plt.plot(np.arange(0, m, 1), rho, label=f't -> {t}')
            plt.xlabel(r"$x$")
            plt.ylabel(r"$\rho$")


            plt.legend()

    # plt.plot(np.arange(0, m, 1), rho)
    # plt.show()
    plt.savefig('diffusion.png')


if __name__ == "__main__":
    sim()
