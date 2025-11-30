"""
Метод Гаусса:

Прямой ход: приведение матрицы к треугольному виду

Обратный ход: последовательное нахождение неизвестных

Метод простой итерации:

Последовательное уточнение решения по формуле: x_i^(k+1) = (b_i - ΣA_ij*x_j^k)/A_ii

Использую значения с предыдущей итерации

Метод Зейделя:

Модификация метода простой итерации

Использую уже вычисленные новые значения на текущей итерации

Более быстрая сходимость

Метод релаксации:

Введение параметра релаксации ω

При ω=1 эквивалентен методу Зейделя

Оптимальный ω ускоряет сходимость
"""

import numpy as np


def gauss_elimination(A, b):
    n = len(A)
    A = A.copy().astype(float)
    b = b.copy().astype(float)

    for i in range(n):
        div = A[i][i]
        A[i] /= div
        b[i] /= div
        for j in range(i + 1, n):
            factor = A[j][i]
            A[j] -= factor * A[i]
            b[j] -= factor * b[i]

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= A[i][j] * x[j]
    return x


def simple_iteration(A, b, x0, E):
    n = len(A)
    k = 0
    x = x0.copy()
    max_error = np.max(np.abs((A @ x - b) / A.diagonal()))

    while max_error >= E:
        x_new = np.zeros(n)
        for i in range(n):
            x_new[i] = (b[i] - np.dot(A[i, :i], x[:i]) -
                        np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
        x = x_new.copy()
        k += 1
        max_error = np.max(np.abs((A @ x - b) / A.diagonal()))

    return x, k, max_error


def seidel_method(A, b, x0, E):
    n = len(A)
    k = 0
    x = x0.copy()
    max_error = np.max(np.abs((A @ x - b) / A.diagonal()))

    while max_error >= E:
        x_new = np.zeros(n)
        for i in range(n):
            x_new[i] = (b[i] - np.dot(A[i, :i], x_new[:i]) -
                        np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
        x = x_new.copy()
        k += 1
        max_error = np.max(np.abs((A @ x - b) / A.diagonal()))

    return x, k, max_error


def relaxation(A, b, E, omega):
    n = len(A)
    x = np.zeros(n)
    x_prev = np.zeros(n)
    k = 0
    r = np.inf

    while r > E:
        for i in range(n):
            x[i] = (1 - omega) * x_prev[i] + (omega / A[i, i]) * \
                   (b[i] - np.sum(A[i, :i] * x[:i]) -
                    np.sum(A[i, i + 1:] * x_prev[i + 1:]))
        r = np.max(np.abs((x - x_prev) / x))
        x_prev = x.copy()
        k += 1

    return x, k, r


if __name__ == "__main__":
    A = np.array([[7, 0.8, 0.9], [0.8, 8, 1], [0.9, 1, 9]])
    b = np.array([63.5, 78.6, 95.3])

    # Метод Гаусса
    x_gauss = gauss_elimination(A, b)
    print("Метод Гаусса:", x_gauss)

    # Метод простой итерации
    x0 = np.zeros(3)
    E = 0.0001
    x_iter, k_iter, error_iter = simple_iteration(A, b, x0, E)
    print("Метод простой итерации:", x_iter)
    print("Количество итераций:", k_iter)
    print("Погрешность:", error_iter)

    # Метод Зейделя
    x_seidel, k_seidel, error_seidel = seidel_method(A, b, x0, E)
    print("Метод Зейделя:", x_seidel)
    print("Количество итераций:", k_seidel)
    print("Погрешность:", error_seidel)

    # Метод релаксации
    print("\nМетод верхней релаксации:")
    for omega in np.arange(0.1, 1.2, 0.1):
        x_relax, k_relax, r_relax = relaxation(A, b, E, omega)
        print(f"Омега = {omega:.1f}: {x_relax}, Итераций: {k_relax}, Погрешность: {r_relax:.6f}")