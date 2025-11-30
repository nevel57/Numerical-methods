"""
Аналитические функции:

U(x, t) - точное решение уравнения

u0(x), u1(t), u2(t) - начальные и граничные условия

f(x, t) - правая часть дифференциального уравнения

Вспомогательные функции:

print_array() - форматированный вывод матриц

print_exact_solution() - вычисление и вывод точного решения

Основной алгоритм solve_for():

Шаг 1: Определение шагов сетки h и τ
Шаг 2: Инициализация расчетной сетки
Шаг 3: Задание граничных и начальных условий
Шаг 4: Решение неявной схемы методом прогонки
Шаг 5: Решение явной схемы
Шаг 6: Вычисление погрешности относительно точного решения
"""


def U(x, t):
    """Точное решение: u(x,t) = x²t² + 3x*t + x + 3"""
    return x * x * t * t + 3 * x * t + x + 3


def u0(x):
    """Начальное условие при t=0: u(x,0) = x + 3"""
    return x + 3


def u1(t):
    """Граничное условие при x=0: u(0,t) = 3"""
    return 3


def u2(t):
    """Граничное условие при x=1: u(1,t) = t² + 3t + 4"""
    return t * t + 3 * t + 4


def f(x, t):
    """Правая часть уравнения: ∂u/∂t - ∂²u/∂x² = f(x,t)"""
    u_t = 2 * x * x * t + 3 * x  # ∂u/∂t
    u_xx = 2 * t * t  # ∂²u/∂x²
    return u_t - u_xx


def print_array(y):
    print("----------")
    for row in y:
        row_str = []
        for k in row:
            row_str.append("%.4f" % k)
        print(" ".join(row_str))
    print("----------")


def print_exact_solution(h, tau, m, n):
    """Печать точного решения"""
    u_f = [[0.0] * (n + 1) for _ in range(m + 1)]
    for j in range(m + 1):
        for i in range(n + 1):
            u_f[j][i] = U(i * h, tau * j)
    print_array(u_f)


def solve_for(sigma, n, m, T):
    h = 1.0 / n
    tau = T * 1.0 / m if sigma != 0 else h * h / 2

    print(f"solving for sigma={sigma} N={n} M={m} T={T}")
    print(f"h={h:.4f} tau={tau:.6f}")

    y = [[0.0] * (n + 1) for _ in range(m + 1)]

    for j in range(m + 1):
        y[j][0] = u1(j * tau)
        y[j][n] = u2(j * tau)
    for i in range(n + 1):
        y[0][i] = u0(i * h)

    z = []

    if sigma != 0:
        for j in range(m):
            alpha = [0.0] * (n + 1)
            beta = [0.0] * (n + 1)

            A = sigma * tau / (h * h)
            B = sigma * tau / (h * h)
            C = 1.0 + 2 * sigma * tau / (h * h)

            for i in range(1, n):
                F = (y[j][i] +
                     (1 - sigma) * tau * (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / (h * h) +
                     tau * f(i * h, j * tau))

                if i == 1:
                    denominator = C
                    alpha[i] = B / denominator
                    beta[i] = (F + A * y[j + 1][0]) / denominator
                else:
                    denominator = C - A * alpha[i - 1]
                    alpha[i] = B / denominator
                    beta[i] = (F + A * beta[i - 1]) / denominator

            for i in range(n - 1, 0, -1):
                if i == n - 1:
                    y[j + 1][i] = alpha[i] * y[j + 1][n] + beta[i]
                else:
                    y[j + 1][i] = alpha[i] * y[j + 1][i + 1] + beta[i]

    else:
        for j in range(m):
            for i in range(1, n):
                phi_ij = f(i * h, j * tau + tau / 2)
                tau_h_h = tau / (h * h)
                y[j + 1][i] = ((1 - 2 * tau_h_h) * y[j][i] +
                               tau_h_h * (y[j][i - 1] + y[j][i + 1]) +
                               tau * phi_ij)

    for j in range(m + 1):
        for i in range(n + 1):
            z.append(abs(U(i * h, j * tau) - y[j][i]))

    print(f"Zmax={max(z):.6f}")
    print("\n")
    return y, max(z)


if __name__ == "__main__":
    print("РЕШЕНИЕ УРАВНЕНИЯ ТЕПЛОПРОВОДНОСТИ")
    print("Точное решение: u(x,t) = x²t² + 3x*t + x + 3")
    print("=" * 60)

    # Тестовые запуски
    print("Точное решение для h=0.1, tau=0.1:")
    print_exact_solution(0.1, 0.1, 10, 10)

    y1, error1 = solve_for(1, 10, 10, 1)
    y2, error2 = solve_for(1, 10, 100, 1)
    y3, error3 = solve_for(0.5, 10, 10, 1)
    y4, error4 = solve_for(0.5, 10, 100, 1)

    print("\n" + "=" * 60)
    print("СРАВНЕНИЕ ПОГРЕШНОСТЕЙ:")
    print(f"sigma=1, N=10, M=10:    {error1:.6f}")
    print(f"sigma=1, N=10, M=100:   {error2:.6f}")
    print(f"sigma=0.5, N=10, M=10:  {error3:.6f}")
    print(f"sigma=0.5, N=10, M=100: {error4:.6f}")
