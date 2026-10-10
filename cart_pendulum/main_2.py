import numpy as np
import matplotlib.pyplot as plt


def calc_acc(
    m: float,
    M: float,
    L: float,
    g: float,
    theta: float,
    theta_1: float,
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


def calc_force(
    theta: float, theta_1: float, omega: float, m: float, M: float, L: float, g: float
) -> float:
    # accelerazione calcolata tramite funzione di guess esponenziale del tipo x(t)=a+(bt+c)e^(-lambda t)
    lam = 1
    theta_2 = -lam * lam * (theta - np.pi / 2) - 2 * lam * theta_1
    x_2 = (m * g * L * np.cos(theta) + 4 * m * L * L * theta_2 / 3) / (
        m * L * np.sin(theta)
    )
    F = (
        (m + M) * x_2
        - m * L * np.sin(theta) * theta_2
        - m * L * theta_1 * theta_1 * np.cos(theta)
    )
    F_max = 100
    if F > F_max:
        F = F_max
    elif F < -F_max:
        F = -F_max
    return F


if __name__ == "__main__":
    dt = 0.001
    omega = 1
    time = 0
    m = 1
    M = 1
    L = 1
    g = 10
    theta = 0.2
    theta_1 = 0
    x = 0
    x_1 = 0
    nSteps = 100000
    theta_vec = np.array([])
    x_vec = np.array([])
    time_vec = np.array([])
    F_vec = np.array([])

    for i in range(nSteps):
        F = calc_force(theta, theta_1, omega, m, M, L, g)
        # F = 0
        acc = calc_acc(m, M, L, g, theta, theta_1, F)
        x_vec = np.append(x_vec, x)
        theta_vec = np.append(theta_vec, theta)
        time_vec = np.append(time_vec, time)
        F_vec = np.append(F_vec, F)
        x_1 += acc[0] * dt
        x += (x_1 + 0.5 * acc[0] * dt) * dt
        theta_1 += acc[1] * dt
        theta += (theta_1 + 0.5 * acc[1] * dt) * dt
        time += dt
    plt.scatter(time_vec, theta_vec)
    # plt.scatter(time_vec, F_vec / 1000)
    plt.show()
    plt.scatter(time_vec, F_vec)
    plt.show()
    plt.scatter(time_vec, x_vec)
    plt.show()
