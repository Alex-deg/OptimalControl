from math import e, cos, sin
import time
import matplotlib.pyplot as plt
from functions import get_valid_input, additional_test
import numpy as np

g = 9.81

# вычисление силы Архимеда (пускай будет)
def compute_archimedes_force(weight):
    return weight * g

# вычисление коэффициента а (пускай будет)
def compute_a(archimedes_force, weight):
    return 1 - archimedes_force / (weight * g)

# функция вывода графиков по map'e 
def print_graphics(x_axis, y_axis, title, x_label, y_label):

    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_axis, 'o-', linewidth=2, markersize=2)
    plt.title(title) # главный заголовок
    plt.xlabel(x_label) # подпись оси х
    plt.ylabel(y_label) # подпись оси у
    plt.tight_layout()
    plt.show()

# основная функция, реализующая численное интегрирование методом Эйлера/Рунге-Кутта 1-го порядка
def solver(h = 0.01, T = 10, control_vector = [1, 1]):

    # задаём начальные условия
    V     = np.zeros(int(T / h))  
    V[0] = 1    
    theta = np.zeros(int(T / h))      
    x     = np.zeros(int(T / h))           
    y     = np.zeros(int(T / h))           

    R, Y = 1200, 2000
    m = 1500
    ro = 1000
    g = 9.81
    Volume = 1.4
    A = ro * g * Volume
    a = 1 - A / (m * g) # принимаем за константу, так как управляющие значения nx, ny, nz

    # получаем значения вектора управления
    nx = control_vector[0]
    ny = control_vector[1]
    
    i = 1

    # цикл численного интегрирования
    while i < T / h:

        V_prev = V[i - 1]
        theta_prev = theta[i - 1]

        # A_x = A * sin(theta_prev)
        # A_y = A * cos(theta_prev)

        # nx = (R - A_x) / (m * g)
        # ny = (Y + A_y) / (m * g)
        
        V[i] = V_prev + h * (a * g * (nx - sin(theta_prev)))
        
        theta[i] = theta_prev + h * (a * g/V_prev * (ny - cos(theta_prev)))
            
        x[i] = x[i - 1] + h * (V_prev * cos(theta_prev)) # угол рыскания = 0 |=> cos(фи) = cos(0) = 1
        y[i] = y[i - 1] + h * (V_prev * sin(theta_prev))
        
        i += 1
    
    return [V, theta, x, y]

if __name__ == "__main__":
    
    # считываем общее количество экспериментов, в которых будем проверять заисимость V, theta, x и y от [nx, ny, nz]
    print('Введите количество экспериментов\nexperiments = ', end='')  
    experiments = input()
    experiments = get_valid_input('int', experiments, [])

    # заголовки для графиков V, theta, x, y соответственно
    titles_for_graphics = [['График изменения скорости V', 'T, с', 'м/с'],
                           ['График изменения угла тангажа', 'T, c', 'радианы'],
                           ['График изменения координаты по оси OX', 'T, c', ''],
                           ['График изменения координаты по оси OY', 'T, c', '']]

    # основной цикл
    for i in range(experiments):

        # ====================ШАГ 1====================

        # вводим вектор управления - веткор нагрузок [nx, ny, nz]
        print(f'Введите значения параметров вектора управления для эксперимента №{i + 1} через запятую')
        print('nx, ny соответственно')
        # единтсвенное нет проверки для вводимых значений [nx, ny, nz]
        control_vector = list(map(float, input().split(', ')))
        # вводим шаг интегрирования
        print('Введите шаг интегрирования\nh = ', end='')
        h = input()
        h = get_valid_input('float', h, [])
        print('Введите время моделирования процесса\nT = ', end='')
        T = input()
        T = get_valid_input('int', T, [])
        # вычисление значений V, theta, x и y спустя некоторое время T, введенное пользователем
        list_of_values = solver(h, T, control_vector) # list_of_values - это V, theta, x, y
        # цикл вывода графиков интересующих нас значений - V, theta, x, y
        x_axis = np.arange(0, T, 0.01)
        for i in range(len(list_of_values)):
            main_title, x_title, y_title = titles_for_graphics[i]
            print_graphics(x_axis, list_of_values[i], main_title, x_title, y_title)
        
        # ====================ШАГ 2====================

        print(f'Введите координаты точки, в которую хотите попасть за {T} секунд:')
        print('x_f, y_f = ', end='')
        
