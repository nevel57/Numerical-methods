"""
Вычисление функции и производной:

Функция: f(x) = Σ x^(k+1) / (2^k * (2k)!)

Аналитическая производная: f'(x) = Σ (k+1) * x^k / (2^k * (2k)!)

Суммирование до точности epsilon = 1e-6

Сравнение методов численного дифференцирования:

Вычисляются три типа разностных производных в точках x = 1.1, 1.2, 1.3, 1.4

Для каждого метода вычисляется погрешность относительно аналитической производной
"""


def f(x):
    return x + (x ** 2) / 4 + (x ** 3) / 96


def right_rectangles(a, b, n, func):
    h = (b - a) / n
    s = 0
    for i in range(1, n + 1):
        s += func(a + i * h)
    return s * h


def trapezoidal(a, b, n, func):
    h = (b - a) / n
    s = (func(a) + func(b)) / 2
    for i in range(1, n):
        s += func(a + i * h)
    return s * h


def simpson(a, b, n, func):
    h = (b - a) / n
    s = func(a) + func(b)
    for i in range(1, n, 2):
        s += 4 * func(a + i * h)
    for i in range(2, n - 1, 2):
        s += 2 * func(a + i * h)
    return s * h / 3


def integrate(a, b, n, func, method, e=0.0001):
    prev_result = method(a, b, n, func)
    n *= 2
    curr_result = method(a, b, n, func)
    while abs(curr_result - prev_result) > e:
        prev_result = curr_result
        n *= 2
        curr_result = method(a, b, n, func)
    return curr_result, n


def true_integral(a, b):
    return (b ** 2) / 2 + (b ** 3) / 12 + (b ** 4) / 96 - ((a ** 2) / 2 + (a ** 3) / 12 + (a ** 4) / 96)


a = 0.2
b = 0.4
n = 10
methods = {
    "прямоугольников": right_rectangles,
    "трапеции": trapezoidal,
    "Симпсона": simpson,
}
F = true_integral(a, b)
for name, method in methods.items():
    result, k = integrate(a, b, n, f, method)
    R = abs(result - F) / abs(F)
    print(f"Метод {name}:")
    print(f"    Приближенное значение интеграла: {result}.6f")
    print(f"    Количество шагов: {k}")
    print(f"    Погрешность: {R}")
