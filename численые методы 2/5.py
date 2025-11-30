"""

"""


def exact_solution(x, t):
    """1. Точное решение для варианта 3: u(x,t) = x³t³ + 3xt + x + 3"""
    return x ** 3 * t ** 3 + 3 * x * t + x + 3


def u0(x):
    """2. Начальное условие: u(x,0) = x + 3"""
    return x + 3


def u1(x):
    """3. Производная начального условия: ∂u/∂t(x,0) = 3x"""
    return 3 * x


def mu1(t):
    """4. Граничное условие при x=0: u(0,t) = 3"""
    return 3


def mu2(t):
    """5. Граничное условие при x=1: u(1,t) = t³ + 3t + 4"""
    return t ** 3 + 3 * t + 4


def f(x, t):
    """6. Правая часть уравнения: f(x,t) = ∂²u/∂t² - ∂²u/∂x² = 6x³t - 6xt³"""
    return 6 * x ** 3 * t - 6 * x * t ** 3


def f_new(x, t, h):
    """7. Модифицированная правая часть с учетом повышенной точности O(h⁴)"""
    f_xxxx = 0

    return f(x, t) - (h ** 2 / 12) * f_xxxx


def solve_wave_equation(sigma, N, M, T):

    h = 1.0 / N
    tau = T / M

    print(f"Решение для sigma={sigma}, N={N}, M={M}, T={T}")
    print(f"h={h:.4f}, tau={tau:.6f}")

    y = [[0.0] * (N + 1) for _ in range(M + 1)]

    for i in range(N + 1):
        y[0][i] = u0(i * h)

    for i in range(1, N):
        y[1][i] = y[0][i] + tau * u1(i * h) + 0.5 * tau ** 2 * (
                (y[0][i + 1] - 2 * y[0][i] + y[0][i - 1]) / h ** 2 +
                f(i * h, 0)
        )

    for j in range(M + 1):
        y[j][0] = mu1(j * tau)
        y[j][N] = mu2(j * tau)

    if sigma == 0:
        for j in range(1, M):
            for i in range(1, N):
                y[j + 1][i] = 2 * y[j][i] - y[j - 1][i] + tau ** 2 * (
                        (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                        f(i * h, j * tau)
                )
    else:
        a = sigma * tau ** 2 / h ** 2
        c = 1 + 2 * a

        for j in range(1, M):
            F = [0.0] * (N + 1)
            for i in range(1, N):
                F[i] = 2 * y[j][i] - y[j - 1][i] + tau ** 2 * (
                        (1 - 2 * sigma) * (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                        sigma * (y[j - 1][i + 1] - 2 * y[j - 1][i] + y[j - 1][i - 1]) / h ** 2 +
                        f(i * h, j * tau)
                )

            F[1] += a * y[j + 1][0]
            F[N - 1] += a * y[j + 1][N]

            alpha = [0.0] * (N + 1)
            beta = [0.0] * (N + 1)

            alpha[1] = a / c
            beta[1] = F[1] / c

            for i in range(2, N):
                denominator = c - a * alpha[i - 1]
                alpha[i] = a / denominator
                beta[i] = (F[i] + a * beta[i - 1]) / denominator

            y[j + 1][N - 1] = (F[N - 1] + a * beta[N - 2]) / (c - a * alpha[N - 2])

            for i in range(N - 2, 0, -1):
                y[j + 1][i] = alpha[i] * y[j + 1][i + 1] + beta[i]

    max_error = 0.0
    for j in range(M + 1):
        for i in range(N + 1):
            exact = exact_solution(i * h, j * tau)
            error = abs(y[j][i] - exact)
            if error > max_error:
                max_error = error

    print(f"Максимальная погрешность: {max_error:.6f}")

    if N == 10 and M == 100:
        print("\nПроверка точности:")
        test_points = [(0, 0), (5, 50), (10, 100), (3, 25), (7, 75)]
        for i, j in test_points:
            x_val = i * h
            t_val = j * tau
            exact_val = exact_solution(x_val, t_val)
            approx_val = y[j][i]
            error_val = abs(exact_val - approx_val)
            print(f"x={x_val:.1f}, t={t_val:.2f}: точн={exact_val:.4f}, прибл={approx_val:.4f}, ошибка={error_val:.4f}")

    print()
    return y, max_error


def solve_wave_equationEX(sigma, N, M, T):
    """Решение с модифицированной правой частью для повышения точности"""
    h = 1.0 / N
    tau = T / M

    print(f"Решение для sigma={sigma}, N={N}, M={M}, T={T}")
    print(f"h={h:.4f}, tau={tau:.6f}")

    y = [[0.0] * (N + 1) for _ in range(M + 1)]

    for i in range(N + 1):
        y[0][i] = u0(i * h)

    for i in range(1, N):
        y[1][i] = y[0][i] + tau * u1(i * h) + 0.5 * tau ** 2 * (
                (y[0][i + 1] - 2 * y[0][i] + y[0][i - 1]) / h ** 2 +
                f_new(i * h, 0, h)
        )

    for j in range(M + 1):
        y[j][0] = mu1(j * tau)
        y[j][N] = mu2(j * tau)

    if sigma == 0:
        for j in range(1, M):
            for i in range(1, N):
                y[j + 1][i] = 2 * y[j][i] - y[j - 1][i] + tau ** 2 * (
                        (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                        f_new(i * h, j * tau, h)
                )
    else:
        a = sigma * tau ** 2 / h ** 2
        c = 1 + 2 * a

        for j in range(1, M):
            F = [0.0] * (N + 1)
            for i in range(1, N):
                F[i] = 2 * y[j][i] - y[j - 1][i] + tau ** 2 * (
                        (1 - 2 * sigma) * (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                        sigma * (y[j - 1][i + 1] - 2 * y[j - 1][i] + y[j - 1][i - 1]) / h ** 2 +
                        f_new(i * h, j * tau, h)
                )

            F[1] += a * y[j + 1][0]
            F[N - 1] += a * y[j + 1][N]

            alpha = [0.0] * (N + 1)
            beta = [0.0] * (N + 1)

            alpha[1] = a / c
            beta[1] = F[1] / c

            for i in range(2, N):
                denominator = c - a * alpha[i - 1]
                alpha[i] = a / denominator
                beta[i] = (F[i] + a * beta[i - 1]) / denominator

            y[j + 1][N - 1] = (F[N - 1] + a * beta[N - 2]) / (c - a * alpha[N - 2])

            for i in range(N - 2, 0, -1):
                y[j + 1][i] = alpha[i] * y[j + 1][i + 1] + beta[i]

    max_error = 0.0
    for j in range(M + 1):
        for i in range(N + 1):
            exact = exact_solution(i * h, j * tau)
            error = abs(y[j][i] - exact)
            if error > max_error:
                max_error = error

    print(f"Максимальная погрешность: {max_error:.6f}")

    if N == 10 and M == 100:
        print("\nПроверка точности:")
        test_points = [(0, 0), (5, 50), (10, 100), (3, 25), (7, 75)]
        for i, j in test_points:
            x_val = i * h
            t_val = j * tau
            exact_val = exact_solution(x_val, t_val)
            approx_val = y[j][i]
            error_val = abs(exact_val - approx_val)
            print(f"x={x_val:.1f}, t={t_val:.2f}: точн={exact_val:.4f}, прибл={approx_val:.4f}, ошибка={error_val:.4f}")

    print()
    return y, max_error


def main():
    """10. ОСНОВНАЯ ПРОГРАММА"""
    print("РЕШЕНИЕ УРАВНЕНИЯ КОЛЕБАНИЙ СТРУНЫ")
    print("Вариант 3: u(x,t) = x³t³ + 3xt + x + 3")
    print("=" * 60)

    print("Явная схема (sigma=0):")
    y1, error1 = solve_wave_equation(0, 10, 10, 1)

    print("Неявная схема (sigma=0.5):")
    y2, error2 = solve_wave_equation(0.5, 10, 10, 1)

    print("Неявная схема (sigma=0.5):")
    y2, error2 = solve_wave_equation(0.5, 100, 100, 1)

    print("Неявная схема (sigma=1):")
    y3, error3 = solve_wave_equation(1, 10, 10, 1)

    # ТЕСТЫ С МОДИФИЦИРОВАННОЙ ПРАВОЙ ЧАСТЬЮ
    print("Новая схема (sigma=0):")
    y1, error1 = solve_wave_equationEX(0, 10, 100, 1)

    print("Новая схема (sigma=0.5):")
    y2, error2 = solve_wave_equationEX(0.5, 10, 100, 1)

    print("Новая схема (sigma=1):")
    y3, error3 = solve_wave_equationEX(1, 10, 100, 1)

    print("ИССЛЕДОВАНИЕ СХОДИМОСТИ:")
    errors = []
    for N, M in [(10, 100), (20, 200), (40, 400)]:
        _, error = solve_wave_equation(0.5, N, M, 1)
        errors.append(error)

    print("\nСходимость при уменьшении шага:")
    for i in range(len(errors) - 1):
        ratio = errors[i] / errors[i + 1]
        print(
            f"h={1 / (10 * 2 ** i):.3f} -> h={1 / (10 * 2 ** (i + 1)):.3f}: погрешность уменьшилась в {ratio:.2f} раз")


if __name__ == "__main__":
    main()