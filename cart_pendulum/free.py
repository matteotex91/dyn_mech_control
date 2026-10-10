import numpy as np
import matplotlib.pyplot as plt


def calc_acc(
    m: float,
    M: float,
    L: float,
    g: float,
    theta: float,
    theta_1: float,
    x: float,
    F: float,
) -> np.ndarray:
    A = np.array(
        [[m + M, -m * L * np.sin(theta)], [-m * L * np.sin(theta), 4 * m * L * L / 3]]
    )
    return np.dot(
        np.linalg.inv(A),
        np.array(
            [m * L * theta_1 * theta_1 * np.cos(theta) + F, -m * g * L * np.cos(theta)]
        ),
    )


if __name__ == "__main__":
    dt = 0.001
    omega = 1
    time = 0
    m = 1
    M = 1
    L = 1
    g = 10
    theta = 0
    theta_1 = 0
    x = 0
    x_1 = 0
    nSteps = 10000
    theta_vec = np.array([])
    x_vec = np.array([])
    time_vec = np.array([])

    for i in range(nSteps):
        F = 0
        acc = calc_acc(m, M, L, g, theta, theta_1, x, F)
        x_vec = np.append(x_vec, x)
        theta_vec = np.append(theta_vec, theta)
        time_vec = np.append(time_vec, time)
        x_1 += acc[0] * dt
        x += (x_1 + 0.5 * acc[0] * dt) * dt
        theta_1 += acc[1] * dt
        theta += (theta_1 + 0.5 * acc[1] * dt) * dt
        time += dt
    plt.scatter(time_vec, theta_vec)
    plt.show()
    plt.scatter(time_vec, x_vec)
    plt.show()
