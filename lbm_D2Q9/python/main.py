import numpy as np
import utils_vortex_shedding as utl
import constants as cst
import matplotlib.pyplot as plt

N_TIMES = 100


def main():
    iscylinder = np.zeros((cst.buffer_ny, cst.buffer_nx), dtype=np.int8)
    utl.generate_isCylinder(iscylinde=iscylinder)

    f = np.zeros((cst.NO_Q, cst.buffer_ny, cst.buffer_nx))
    feq = f.copy()
    f[:] = feq[:]

    rho = np.ones((cst.buffer_ny, cst.buffer_nx), dtype=np.float64)

    ux = np.zeros((cst.buffer_ny, cst.buffer_nx), dtype=np.float64)
    uy = np.zeros((cst.buffer_ny, cst.buffer_nx), dtype=np.float64)

    for time in range(N_TIMES):
        utl.calculate_f_equilibrium(
            feq=feq, iscylinde=iscylinder, rho=rho, ux=ux, uy=uy
        )

        utl.collision(f=f, feq=feq, iscylinde=iscylinder)

        utl.apply_zou_he_boundary(f=f, rho=rho)
        utl.apply_wall_bounceback(f=f, rho=rho)
        utl.apply_obstacle_bounceback(f=f, iscylinder=iscylinder, rho=rho)

        f = utl.streaming(f=f, iscylinde=iscylinder)

        utl.calculate_microsopic_quantities(
            f=f, iscylinde=iscylinder, rho=rho, ux=ux, uy=uy
        )
    np.savetxt('data.out',rho,delimiter=' ')
    plt.imshow(iscylinder)

    plt.imshow(rho)
    plt.show()


if __name__ == "__main__":
    # plt.imshow(iscylinder)
    # plt.colorbar()
    # plt.show()
    main()
