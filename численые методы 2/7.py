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


def thomas_algorithm(a, b, c, d, n):
    """Метод прогонки для трехдиагональной системы"""
    # Прямой ход
    c_prime = [0.0] * n
    d_prime = [0.0] * n

    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * c_prime[i - 1]
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom

    # Обратный ход
    x = [0.0] * n
    x[n - 1] = d_prime[n - 1]

    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x


def solve_wave_equation(sigma, N, M, T):
    """8. Решение уравнения колебаний струны σ-методом"""
    h = 1.0 / N
    tau = T / M

    print(f"Решение для sigma={sigma}, N={N}, M={M}, T={T}")
    print(f"h={h:.4f}, tau={tau:.6f}")

    # Инициализация сетки
    y = [[0.0] * (N + 1) for _ in range(M + 1)]

    # Начальные условия
    for i in range(N + 1):
        y[0][i] = u0(i * h)

    # Первый временной слой (используем явную схему)
    for i in range(1, N):
        y[1][i] = y[0][i] + tau * u1(i * h) + 0.5 * tau ** 2 * (
                (y[0][i + 1] - 2 * y[0][i] + y[0][i - 1]) / h ** 2 +
                f(i * h, 0)
        )

    # Граничные условия
    for j in range(M + 1):
        y[j][0] = mu1(j * tau)
        y[j][N] = mu2(j * tau)

    # Основной цикл по времени
    if sigma == 0:
        # ЯВНАЯ СХЕМА
        for j in range(1, M):
            for i in range(1, N):
                y[j + 1][i] = 2 * y[j][i] - y[j - 1][i] + tau ** 2 * (
                        (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                        f(i * h, j * tau)
                )
    else:
        # НЕЯВНАЯ СХЕМА - ИСПРАВЛЕННАЯ ВЕРСИЯ
        a_coef = sigma * tau ** 2 / h ** 2
        b_coef = 1 + 2 * a_coef

        for j in range(1, M):
            # Подготовка данных для метода прогонки
            n_inner = N - 1  # количество внутренних узлов
            a = [a_coef] * n_inner  # нижняя диагональ
            b = [b_coef] * n_inner  # главная диагональ
            c = [a_coef] * n_inner  # верхняя диагональ
            d = [0.0] * n_inner  # правая часть

            # Заполнение правой части
            for i in range(1, N):
                k = i - 1  # индекс в массивах прогонки
                d[k] = (2 * y[j][i] - y[j - 1][i] +
                        tau ** 2 * (
                                (1 - sigma) * (y[j][i + 1] - 2 * y[j][i] + y[j][i - 1]) / h ** 2 +
                                f(i * h, j * tau)
                        ))

            # Учет граничных условий в правой части
            d[0] += a_coef * y[j + 1][0]  # левая граница
            d[-1] += a_coef * y[j + 1][N]  # правая граница

            # Решение системы методом прогонки
            y_inner = thomas_algorithm(a, b, c, d, n_inner)

            # Запись результатов
            for i in range(1, N):
                y[j + 1][i] = y_inner[i - 1]

    # Вычисление погрешности
    max_error = 0.0
    for j in range(M + 1):
        for i in range(N + 1):
            exact = exact_solution(i * h, j * tau)
            error = abs(y[j][i] - exact)
            if error > max_error:
                max_error = error

    print(f"Максимальная погрешность: {max_error:.6f}")

    # Проверка точности в характерных точках
    if M >= 100:
        print("\nПроверка точности:")
        test_points = [(0, 0), (N // 2, M // 2), (N, M), (N // 4, M // 4), (3 * N // 4, 3 * M // 4)]
        for i, j in test_points:
            if j < len(y) and i < len(y[0]):
                x_val = i * h
                t_val = j * tau
                exact_val = exact_solution(x_val, t_val)
                approx_val = y[j][i]
                error_val = abs(exact_val - approx_val)
                print(
                    f"x={x_val:.1f}, t={t_val:.2f}: точн={exact_val:.4f}, прибл={approx_val:.4f}, ошибка={error_val:.6f}")

    print()
    return y, max_error


def main():
    """ОСНОВНАЯ ПРОГРАММА"""
    print("РЕШЕНИЕ УРАВНЕНИЯ КОЛЕБАНИЙ СТРУНЫ")
    print("Вариант 3: u(x,t) = x³t³ + 3xt + x + 3")
    print("=" * 60)

    # ТЕСТИРОВАНИЕ РАЗНЫХ СХЕМ
    print("Явная схема (sigma=0) - грубая сетка:")
    y1, error1 = solve_wave_equation(0, 10, 10, 1)

    print("Явная схема (sigma=0) - мелкая сетка:")
    y2, error2 = solve_wave_equation(0, 10, 100, 1)

    print("Схема Кранка-Николсон (sigma=0.5) - грубая сетка:")
    y3, error3 = solve_wave_equation(0.5, 10, 10, 1)

    print("Схема Кранка-Николсон (sigma=0.5) - мелкая сетка:")
    y4, error4 = solve_wave_equation(0.5, 10, 100, 1)

    print("Неявная схема (sigma=1) - грубая сетка:")
    y5, error5 = solve_wave_equation(1, 10, 10, 1)

    print("Неявная схема (sigma=1) - мелкая сетка:")
    y6, error6 = solve_wave_equation(1, 10, 100, 1)

    # АНАЛИЗ РЕЗУЛЬТАТОВ
    print("\n" + "=" * 60)
    print("СРАВНЕНИЕ РЕЗУЛЬТАТОВ:")
    print("=" * 60)
    print(f"Явная схема, M=10:    {error1:.6f}")
    print(f"Явная схема, M=100:   {error2:.6f}")
    print(f"Схема Кранка-Николсон, M=10:  {error3:.6f}")
    print(f"Схема Кранка-Николсон, M=100: {error4:.6f}")
    print(f"Неявная схема, M=10:  {error5:.6f}")
    print(f"Неявная схема, M=100: {error6:.6f}")

    # АНАЛИЗ СХОДИМОСТИ
    print("\n" + "=" * 60)
    print("АНАЛИЗ СХОДИМОСТИ ДЛЯ ЯВНОЙ СХЕМЫ:")
    print("=" * 60)

    errors_explicit = []
    for N, M in [(10, 100), (20, 200), (40, 400)]:
        _, error = solve_wave_equation(0, N, M, 1)
        errors_explicit.append(error)
        print(f"N={N}, M={M}: погрешность = {error:.6f}")

    print("\nПорядок сходимости явной схемы:")
    for i in range(len(errors_explicit) - 1):
        ratio = errors_explicit[i] / errors_explicit[i + 1]
        order = 1.0  # Ожидаемый порядок для волнового уравнения
        print(f"h уменьшился в 2 раза: погрешность уменьшилась в {ratio:.2f} раз")


if __name__ == "__main__":
    main()