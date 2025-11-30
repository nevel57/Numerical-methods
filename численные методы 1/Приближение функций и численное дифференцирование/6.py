"""
Вычисление функции:

Используется разложение в ряд: f(x) = Σ x^(k+1) / (2^k * (2k)!)

Суммирование до достижения точности epsilon = 1e-6

Численное дифференцирование:

Вычисляю производные в точках x = 1.025, 1.05, 1.075, 1.10

Используются три метода:

Правая разность: f'(x) ≈ [f(x+h) - f(x)] / h

Левая разность: f'(x) ≈ [f(x) - f(x-h)] / h

Центральная разность: f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
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


a = 1
h = 0.1
i_values = [0.25 * i for i in range(1, 5)]  # i = 0.25, 0.5, 0.75, 1.0
x_values = [a + i * h for i in i_values]

print("Таблица разностных производных:")
print("-" * 80)
print(" x_i   |   f(x_i)   |   f_x (правая) |   f_y (левая) |   f_0 (центральная)")
print("-" * 80)

for x in x_values:
    f_x = (compute_f(x + h) - compute_f(x)) / h  # Правая разность
    f_y = (compute_f(x) - compute_f(x - h)) / h  # Левая разность
    f_0 = (compute_f(x + h) - compute_f(x - h)) / (2 * h)  # Центральная разность

    print(f" {x:.2f}  |  {compute_f(x):.6f}  |  {f_x:.6f}      |  {f_y:.6f}     |  {f_0:.6f}")

print("-" * 80)
