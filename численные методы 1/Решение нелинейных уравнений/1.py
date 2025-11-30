"""
Метод простой итерации:

Преобразование уравнения к виду x = g(x)

Итерационная формула: xₙ₊₁ = 2 - sin(1/xₙ)

Условие остановки по относительной погрешности

Метод Ньютона:

Использование производной для ускорения сходимости

Итерационная формула: xₙ₊₁ = xₙ - f(xₙ)/f'(xₙ)

Более быстрая квадратичная сходимость

Метод деления отрезка пополам:

Требует начальный интервал [a,b] с разными знаками f(a) и f(b)

Последовательное сужение интервала в 2 раза

Гарантированная сходимость
"""

import math


def f(x):
    return x - 2 + math.sin(1 / x)


def simple_iteration_method(x0, epsilon):
    k = 0
    x_prev = x0

    while True:
        x_next = 2 - math.sin(1 / x_prev)
        k += 1

        if abs((x_next - x_prev) / x_next) < epsilon:
            break

        print(f"Итерация {k}: x = {x_next:.6f}, невязка = {abs(f(x_next)):.6f}")
        x_prev = x_next

    return x_next, k


def newton_method(x0, epsilon):
    def df(x):
        return 1 - math.cos(1 / x) / (x ** 2)

    k = 0
    x_prev = x0
    x_next = x_prev - f(x_prev) / df(x_prev)

    while abs((x_next - x_prev) / x_next) >= epsilon:
        k += 1
        x_prev = x_next
        x_next = x_prev - f(x_prev) / df(x_prev)
        print(f"Итерация {k}: x = {x_next:.6f}, невязка = {abs(f(x_next)):.6f}")

    return x_next, k


def bisection_method(a, b, epsilon, max_iter=100):
    if f(a) * f(b) >= 0:
        raise ValueError("Функция должна иметь разные знаки на концах интервала")

    k = 0
    x_prev = a
    x_next = (a + b) / 2

    while abs((x_next - x_prev) / x_next) > epsilon and k < max_iter:
        k += 1
        x_prev = x_next

        if f(a) * f(x_next) < 0:
            b = x_next
        else:
            a = x_next

        x_next = (a + b) / 2
        print(f"Итерация {k}: x = {x_next:.6f}, невязка = {abs(f(x_next)):.6f}")

    return x_next, k


if __name__ == "__main__":
    epsilon = 0.0001

    print("А) Метод простой итерации:")
    print("=" * 50)
    root1, iter1 = simple_iteration_method(1.2, epsilon)
    print(f"Результат: x = {root1:.6f}, итераций: {iter1}")
    print(f"Невязка: {abs(f(root1)):.6f}\n")

    print("Б) Метод Ньютона:")
    print("=" * 50)
    root2, iter2 = newton_method(1.0, epsilon)
    print(f"Результат: x = {root2:.6f}, итераций: {iter2}")
    print(f"Невязка: {abs(f(root2)):.6f}\n")

    print("В) Метод деления отрезка пополам:")
    print("=" * 50)
    try:
        root3, iter3 = bisection_method(1.0, 2.0, epsilon)
        print(f"Результат: x = {root3:.6f}, итераций: {iter3}")
        print(f"Невязка: {abs(f(root3)):.6f}")
    except ValueError as e:
        print(f"Ошибка: {e}")