import numpy as np
import matplotlib.pyplot as plt

x = 0
a = 0
t = 0

x_target = 1
a_target = 0

x1 = 0
a1 = 0

dt = 1e-3


nsteps = 5000
tfin = nsteps * dt

const_m = 1
const_M = 1
const_L = 1
const_g = 9.806

fourier_N = 5
fourier_F0 = 0.01
fourier_a = np.array([0, 0, 0, 0, 0])
fourier_b = np.array([0, 0, 0, 0, 0])
fourier_DELTA = 0.1


def force(t: float) -> float:
    fsum = fourier_F0
    for i in range(fourier_N):
        omegan = np.pi * (i + 1) / tfin
        fsum += fourier_a[i] * np.cos(omegan * t) + fourier_b[i] * np.sin(omegan * t)
    return fsum


def acc_matrix(a: float) -> np.ndarray:
    return np.linalg.inv(
        np.array(
            [
                [const_M + const_m, -0.5 * const_m * const_L * np.cos(a)],
                [-0.5 * const_m * const_L * np.cos(a), const_m * const_L * const_L / 3],
            ]
        )
    )


def acc_vector(a: float, a1: float, t: float) -> np.ndarray:
    return np.array(
        [
            0.5 * const_m * const_L * a1 * a1 * np.sin(a) + force(t),
            -0.5 * const_m * const_L * np.sin(a) * const_g,
        ]
    )


# x_vec = []
# a_vec = []
# x1_vec = []
# a1_vec = []
# t_vec = []

# for i in range(nsteps):
#     x_vec.append(x)
#     a_vec.append(a)
#     t_vec.append(t)
#     x1_vec.append(x1)
#     a1_vec.append(a1)
#     acc_mat = acc_matrix(a)
#     acc_vec = acc_vector(a, a1, t)
#     acc = np.dot(acc_mat, acc_vec)
#     x1 += acc[0] * dt
#     a1 += acc[1] * dt
#     x += (0.5 * acc[0] * dt + x1) * dt
#     a += (0.5 * acc[1] * dt + a1) * dt
#     t += dt


# plt.plot(t_vec, x_vec)
# plt.plot(t_vec, a_vec)
# plt.show()


def err_fun(t0, x0, a0, tf, xf, af, nsteps) -> float:
    x = x0
    x1 = 0
    a = a0
    a1 = 0
    t = t0
    dt = (tf - t0) / nsteps
    for i in range(nsteps):
        acc_mat = acc_matrix(a)
        acc_vec = acc_vector(a, a1, t)
        acc = np.dot(acc_mat, acc_vec)
        x1 += acc[0] * dt
        a1 += acc[1] * dt
        x += (0.5 * acc[0] * dt + x1) * dt
        a += (0.5 * acc[1] * dt + a1) * dt
        t += dt
    return np.sqrt(
        np.square(x - xf) + np.square(a - af) + np.square(a1) + np.square(x1)
    )


error = err_fun(0, 0, 0, 1, 1, 0, 1000)
print(error)
