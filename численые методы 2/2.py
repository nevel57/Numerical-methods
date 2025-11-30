"""
Постановка задачи:

Дифференциальное уравнение: y'' - 3x·y = -f(x)

Граничные условия: y(0) = 3, y(1) = 13/3

Аналитическое решение: y(x) = x³ + x/3 + 3

Метод конечных разностей:

Разбиение интервала [0,1] на N равных частей

Аппроксимация второй производной разностной схемой

Построение трехдиагональной системы уравнений

Метод прогонки:

Прямой ход: вычисление прогоночных коэффициентов α и β

Обратный ход: последовательное нахождение решения
"""


def f(x):
    return 3 * x ** 4 + x ** 2 + 3 * x


N = int(input("Введите N: "))
h = 1 / N
n1 = 0
n2 = 7
A = 1 / h ** 2
B = A

x = [0.0] * (N + 1)
y = [0.0] * (N + 1)
alpha = [0.0] * (N + 1)
beta = [0.0] * (N + 1)

for i in range(0, N + 1):
    x[i] = i * h

for i in range(0, N + 1):
    q = 3 * x[i]
    C = -2 / h ** 2 - q
    F = -f(x[i])

    if i == 0:
        alpha[0] = 0
        beta[0] = 3.0
    elif i == N:
        pass
    else:
        denominator = C + A * alpha[i - 1]
        alpha[i] = -B / denominator
        beta[i] = (F - A * beta[i - 1]) / denominator

y[N] = 13.0 / 3
for i in range(N - 1, -1, -1):
    y[i] = alpha[i] * y[i + 1] + beta[i]


def exact_solution(x):
    return x ** 3 + x / 3 + 3


print("\nРЕЗУЛЬТАТЫ:")
print("x\t\tY_прибл\t\tU_точн\t\tПогрешность")
print("-" * 50)

max_error = 0.0
for i in range(N + 1):
    U_exact = exact_solution(x[i])
    error = abs(y[i] - U_exact)
    if error > max_error:
        max_error = error

    if N <= 10 or i % (N // 10) == 0 or i == N:
        print(f"{x[i]:.4f}\t\t{y[i]:.8f}\t\t{U_exact:.8f}\t\t{error:.8f}")

print(f"\nМаксимальная погрешность: {max_error:.16f}")

if N <= 20:
    print("\nПОЛНАЯ ТАБЛИЦА:")
    print("i\tx\t\tY_прибл\t\tU_точн\t\tПогрешность")
    print("-" * 70)
    for i in range(N + 1):
        U_exact = exact_solution(x[i])
        error = abs(y[i] - U_exact)
        print(f"{i}\t{x[i]:.4f}\t\t{y[i]:.8f}\t\t{U_exact:.8f}\t\t{error:.8f}")