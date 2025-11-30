

def exact_solution(x):
    """Точное решение: U(x) = 3x² - 2x"""
    return 3 * x * x - 2 * x


def f(x, U):
    """Правая часть уравнения: U' = (3x² + U)/x"""
    return (3 * x * x + U) / x


def euler_step(x, y, h):
    """Метод Эйлера: y_{n+1} = y_n + h*f(x_n, y_n)"""
    return y + h * f(x, y)


def rk2_step(x, y, h):
    """Метод Рунге-Кутта 2-го порядка (модифицированный Эйлер)"""
    k1 = f(x, y)
    k2 = f(x + h, y + h * k1)
    return y + (h / 2) * (k1 + k2)


def rk4_step(x, y, h):
    """Метод Рунге-Кутта 4-го порядка"""
    k1 = f(x, y)
    k2 = f(x + h / 2, y + (h / 2) * k1)
    k3 = f(x + h / 2, y + (h / 2) * k2)
    k4 = f(x + h, y + h * k3)
    return y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def solve_ode(method_name, method_func, h, x_end=2.0):
    print(f"{method_name}")
    print(f"h = {h}")
    print("x\t\tY_прибл\t\tU_точн\t\tПогрешность")
    print("-" * 50)

    x = 1.0
    y = 1.0
    max_error = 0.0

    # Выводим начальную точку
    U_exact = exact_solution(x)
    error = abs(y - U_exact)
    max_error = error
    print(f"{x:.2f}\t\t{y:.6f}\t\t{U_exact:.6f}\t\t{error:.6f}")

    # Количество шагов
    N = int((x_end - 1.0) / h)

    for i in range(N):
        y = method_func(x, y, h)
        x += h

        U_exact = exact_solution(x)
        error = abs(y - U_exact)
        max_error = max(max_error, error)

        print(f"{x:.2f}\t\t{y:.6f}\t\t{U_exact:.6f}\t\t{error:.6f}")

    print(f"Максимальная погрешность: {max_error:.6f}\n")
    return max_error


def test_exact_solution():
    """Тестируем точное решение"""
    print("Проверка точного решения:")
    print("x\tU(x) = 3x² - 2x\tU'(x)\t(3x²+U)/x")
    for x in [1.0, 1.5, 2.0]:
        U = exact_solution(x)
        U_derivative_exact = 6 * x - 2  # Производная от 3x² - 2x
        U_derivative_eq = (3 * x * x + U) / x  # Из уравнения
        print(f"{x}\t{U:.4f}\t\t{U_derivative_exact:.4f}\t{U_derivative_eq:.4f}")


def main():
    print("РЕШЕНИЕ ЗАДАЧИ КОШИ")
    print("Уравнение: U' = (3x² + U)/x")
    print("Начальное условие: U(1) = 1")
    print("Точное решение: U(x) = 3x² - 2x")
    print("=" * 60)

    # Проверяем точное решение
    test_exact_solution()
    print("\n" + "=" * 60)

    METHODS = {
        "МЕТОД ЭЙЛЕРА": euler_step,
        "МЕТОД РУНГЕ-КУТТА 2-го ПОРЯДКА": rk2_step,
        "МЕТОД РУНГЕ-КУТТА 4-го ПОРЯДКА": rk4_step
    }

    h_values = [0.1, 0.05]
    results = {}

    for h in h_values:
        print(f"\nШАГ h = {h}")
        print("=" * 50)

        h_results = {}
        for method_name, method_func in METHODS.items():
            error = solve_ode(method_name, method_func, h)
            h_results[method_name] = error

        results[h] = h_results

    # Вывод сравнения
    print("\n" + "=" * 70)
    print("СРАВНЕНИЕ ПОГРЕШНОСТЕЙ")
    print("=" * 70)
    print("Метод\t\t\t|z|_c (h=0.1)\t|z|_c (h=0.05)\tУменьшение")
    print("-" * 70)

    for method_name in METHODS.keys():
        error_h1 = results[0.1][method_name]
        error_h2 = results[0.05][method_name]
        reduction = error_h1 / error_h2 if error_h2 > 1e-10 else float('inf')

        short_name = method_name.replace("МЕТОД ", "")
        print(f"{short_name}\t\t{error_h1:.6f}\t{error_h2:.6f}\t{reduction:.2f} раз")


if __name__ == "__main__":
    main()