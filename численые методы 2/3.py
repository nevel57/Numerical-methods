"""
Постановка задачи:

Уравнение: -Δu = f(x₁,x₂) в области [0,1]×[0,1]

Точное решение: u(x₁,x₂) = 3x₁³ + x₂³ + 3x₁ + x₂ + 3

Правая часть: f(x₁,x₂) = -18x₁ - 6x₂

Разностная схема:

Аппроксимация лапласиана пятиточечным шаблоном

Итерационный метод верхней релаксации

Параметр релаксации ω для ускорения сходимости

Метод решения:

Построение равномерной сетки N×N

Задание граничных условий из точного решения

Итерационное уточнение решения во внутренних узлах
"""


def solve_poisson_equation():
    """Решение уравнения Пуассона методом верхней релаксации"""

    N = int(input("N: "))
    epsilon = float(input("epsilon: "))

    def u_exact(x1, x2):
        """Точное решение уравнения Пуассона"""
        return 3 * x1 ** 3 + x2 ** 3 + 3 * x1 + x2 + 3

    def f(x1, x2):
        """Правая часть уравнения Пуассона: -Δu = -f(x)"""
        return -18 * x1 - 6 * x2

    def solve(omega):
        """Решение уравнения Пуассона методом верхней релаксации"""
        h = 1.0 / N  # Шаг сетки
        u = [[0.0] * (N + 1) for _ in range(N + 1)]

        for i in range(N + 1):
            for j in range(N + 1):
                if i == 0 or i == N or j == 0 or j == N:
                    u[i][j] = u_exact(i * h, j * h)

        k = 0
        while True:
            max_diff = 0.0

            for i in range(1, N):
                for j in range(1, N):
                    old = u[i][j]
                    new = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] + h * h * f(i * h, j * h)) / 4
                    u[i][j] = omega * new + (1 - omega) * old
                    max_diff = max(max_diff, abs(u[i][j] - old))

            k += 1
            if max_diff < epsilon or k > 10000:
                break

        max_err = 0.0
        for i in range(1, N):
            for j in range(1, N):
                err = abs(u[i][j] - u_exact(i * h, j * h))
                max_err = max(max_err, err)

        return k, max_err, u

    print("\nОПТИМАЛЬНЫЙ ПАРАМЕТР РЕЛАКСАЦИИ:")
    print("ω\tИтерации\tПогрешность")
    print("-" * 30)

    results = []
    for omega in [1.0, 1.2, 1.5, 1.7, 1.9]:
        k, err, _ = solve(omega)
        results.append((omega, k, err))
        print(f"{omega}\t{k}\t\t{err:.2e}")

    best_omega, best_iter, best_error = min(results, key=lambda x: x[1])
    print(f"\nω_opt = {best_omega} (минимум итераций: {best_iter})")

    print("\nРЕЗУЛЬТАТЫ ДЛЯ ω = 1.7 (метод Зейделя):")
    iterations, error, u = solve(1.7)
    print(f"Количество итераций: {iterations}")
    print(f"Максимальная погрешность: {error:.2e}")

    print("\nПроверка в точках:")
    print("i\tj\tx1\tx2\tY_прибл\t\tU_точн\t\tПогрешность")
    print("-" * 65)

    for i in [N // 4, N // 2, 3 * N // 4]:
        for j in [N // 4, N // 2, 3 * N // 4]:
            x1 = i * (1.0 / N)
            x2 = j * (1.0 / N)
            exact = u_exact(x1, x2)
            approx = u[i][j]
            error = abs(approx - exact)
            print(f"{i}\t{j}\t{x1:.2f}\t{x2:.2f}\t{approx:.6f}\t{exact:.6f}\t{error:.2e}")


if __name__ == "__main__":
    solve_poisson_equation()