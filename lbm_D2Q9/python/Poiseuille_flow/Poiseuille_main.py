import Poiseuille_flow as utl
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import constants as cst
import datetime

obj = datetime.datetime.now()
obj.strftime("%B_%d_%Y_%I_%M_%p")


# initialization
#


def main():
    ntime = 10000

    rho = np.ones((utl.NY, utl.NX), dtype=np.float64)
    ux = np.zeros((utl.NY, utl.NX), dtype=np.float64)
    uy = np.zeros((utl.NY, utl.NX), dtype=np.float64)

    feq = np.zeros((cst.NO_Q, utl.NY, utl.NX), dtype=np.float64)
    f = feq.copy()

    utl.calculate_f_eq(feq=feq, rho=rho, ux=ux, uy=uy)
    for k in range(cst.NO_Q):
        f[k] = cst.w[k] * rho

    # f = feq.copy()

    for t in tqdm(range(ntime)):
        utl.calculate_macroscopic_quan(f=f, rho=rho, ux=ux, uy=uy)

        utl.calculate_f_eq(feq=feq, rho=rho, ux=ux, uy=uy)

        utl.collision(f=f, feq=feq, ux=ux, uy=uy)
        f = utl.streaming(f=f)
        # if t % 100 == 0:

        #     print(f"{t}/n", feq)

    plt.plot(np.arange(0, utl.NY, 1), uy[:, utl.NX//2])
    plt.show()

    u = ux ** 2 + uy ** 2
    u_mag = u ** 0.5
    fig, ax = plt.subplots()
    im = ax.imshow(u_mag)
    # plt.colorbar(im)

    ax.quiver(ux.T, uy.T)
    # print(uy)

    obj = datetime.datetime.now()
    filename = obj.strftime("%B_%d_%Y_%I_%M_%p")
    fig.savefig(f'./results//images//{filename}.png')


if __name__ == "__main__":
    main()
