"""
Задаю точку вычисления x = 0.1 и точность e = 0.0001

Использую разложение в степенной ряд для производной

Вычисляю каждый член ряда по формуле: (k+1) * x^k / (2^k * (2k)!)

Суммирую члены, пока абсолютное значение очередного члена не станет меньше заданной точности

Возвращаю полученную сумму как значение производной
"""

import math


def compute_derivative(x: float, e: float) -> float:
    sum_derivative = 0.0
    k = 0
    while True:
        term = (k + 1) * (x ** k) / ((2 ** k) * math.factorial(2 * k))

        sum_derivative += term

        if abs(term) < e:
            break

        k += 1

    return sum_derivative


x = 0.1
e = 0.0001

derivative = compute_derivative(x, e)

print(f"f'({x}) ≈ {derivative:.6f} (с точностью {e})")