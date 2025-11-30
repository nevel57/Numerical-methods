"""
Вычисление функции и производной:

Функция: f(x) = Σ x^(k+1) / (2^k * (2k)!)

Аналитическая производная: f'(x) = Σ (k+1) * x^k / (2^k * (2k)!)

Суммирование до точности epsilon = 1e-6

Сравнение методов численного дифференцирования:

Вычисляю три типа разностных производных в точках x = 1.1, 1.2, 1.3, 1.4

Для каждого метода вычисляется погрешность относительно аналитической производной
"""

import math


def compute_f(x, epsilon=1e-6):
    sum_f = 0.0
    k = 0
    while True:
        term = (x ** (k + 1)) / ((2 ** k) * math.factorial(2 * k))
        sum_f += term
        if abs(term) < epsilon:
            break
        k += 1
    return sum_f


def compute_derivative(x, epsilon=1e-6):
    sum_deriv = 0.0
    k = 0
    while True:
        term = (k + 1) * (x ** k) / ((2 ** k) * math.factorial(2 * k))
        sum_deriv += term
        if abs(term) < epsilon:
            break
        k += 1
    return sum_deriv


h = 0.1
a = 1
points = [a + i * h for i in range(1, 5)]  # Точки x_i = 0.2, 0.3, 0.4, 0.5
print("Погрешности разностных производных:")
print("-" * 130)
print(
    "  x_i   |  f'(x_i)  |  Правая F(x) |  Погрешность z(x) |"
    "  Левая F_x |  Погрешность z(x) |  Центральная Fx0 |  Погрешность z(x)"
)
print("-" * 130)
for x in points:
    f_prime = compute_derivative(x)
    F_right = (compute_f(x + h) - compute_f(x)) / h
    z_right = abs(f_prime - F_right)
    F_left = (compute_f(x) - compute_f(x - h)) / h
    z_left = abs(f_prime - F_left)
    F_central = (compute_f(x + h) - compute_f(x - h)) / (2 * h)
    z_central = abs(f_prime - F_central)

    print(
        f" {x:.1f}    |  {f_prime:.6f} |  {F_right:.6f}    |  {z_right:.6f}         |  "
        f"{F_left:.6f}  |  {z_left:.6f}         |  {F_central:.6f}        |  {z_central:.6f}")

print("-" * 130)
