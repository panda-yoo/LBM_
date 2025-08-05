import numpy as np
import matplotlib.pyplot as plt


def show_data():
    data = np.genfromtxt('cmake-build-debug\\data.dat', delimiter=' ', skip_header=1)
    data = data.reshape(100, 51, 3)
    t = 3

    # t = data[0,:,]
    x = data[t, :, 1]
    rho = data[t, :, 2]

    x = data[0, :, 1]
    rho = data[0, :, 2]
    plt.plot(x, rho, c='r', label="at t = 0")

    x = data[-1, :, 1]
    rho = data[-1, :, 2]
    plt.plot(x, rho, c='b', label=f"at t = {data.shape[0]}")
    plt.legend()
    plt.savefig("data//d1q3.png")
    plt.show()


def show_data1():
    data2 = np.genfromtxt('cmake-build-debug\\data_d2q9.dat', delimiter=' ', skip_header=1)
    x = data2[:, 0]
    y = data2[:, 1]
    rho = data2[:, 2]

    x, y = np.meshgrid(x, y)

    n = (int)((data2.shape[0]) ** 0.5)

    rho = rho.reshape(n, n)
    print(data2.shape)

    im = plt.imshow(rho)
    plt.colorbar(mappable=im)
    plt.show()


if __name__ == '__main__':
    show_data()
    pass
