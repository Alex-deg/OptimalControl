from math import e, cos, sin
import time
import matplotlib.pyplot as plt
from functions import get_valid_input, additional_test
import numpy as np
from scipy.optimize import minimize

g = 9.81

# вычисление силы Архимеда
def compute_archimedes_force(weight):
    return weight * g

# вычисление коэффициента а
def compute_a(archimedes_force, weight):
    return 1 - archimedes_force / (weight * g)

# функция вывода графиков
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

        V[i] = V_prev + h * (a * g * (nx - sin(theta_prev)))    
        theta[i] = theta_prev + h * (a * g/V_prev * (ny - cos(theta_prev)))
        x[i] = x[i - 1] + h * (V_prev * cos(theta_prev)) # угол рыскания = 0 |=> cos(фи) = cos(0) = 1
        y[i] = y[i - 1] + h * (V_prev * sin(theta_prev))
        
        i += 1
    
    return [V, theta, x, y]

# first step
def fixed_control():
    # вводим вектор управления - веткор нагрузок [nx, ny, nz]
    print(f'Введите значения параметров вектора управления для эксперимента через запятую')
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
    print_graphics(list_of_values[2], list_of_values[3], 'График зависимости y от x', 'x', 'y')

# second step
def optimal_control():
    print(f'Введите координаты точки, в которую хотите попасть:')
    print('x_f, y_f = ', end='')
    x_f, y_f = map(int, input().split(', '))
    print(f'Введите время за которое хотите добраться из [0; 0] в [{x_f}; {y_f}]')
    T = input()
    T = get_valid_input('int', T, [])
    p0_optimal, success = simplex_optimization(x_f, y_f, T)
    if success:
        print(f"Оптимальные начальные условия:")
        print(f"  p_V(0)   = {p0_optimal[0]:.6f}")
        print(f"  p_theta(0) = {p0_optimal[1]:.6f}")
        print(f"  p_x(0)   = {p0_optimal[2]:.6f}")
        print(f"  p_y(0)   = {p0_optimal[3]:.6f}")
        
        # 2. Интегрируем систему с оптимальными параметрами
        time, state, controls = integrate_optimal_system(p0_optimal, T)
        
        # 3. Анализируем результаты
        final_state = state[-1]
        print(f"\nКонечное состояние:")
        print(f"  x(T) = {final_state[2]:.2f} м (цель: {x_f} м)")
        print(f"  y(T) = {final_state[3]:.2f} м (цель: {y_f} м)")
        print(f"  p_V(T) = {final_state[4]:.6f} (должно быть 0)")
        print(f"  p_theta(T) = {final_state[5]:.6f} (должно быть 0)")

        # 4. Отрисовка графиков
        plot_optimal_results_simple(time, state, controls, x_f, y_f)
    else:
        print('Оптимизация не сошлась(')


def integrate_optimal_system(p0, T, dt=0.05):
 
    # Параметры системы
    g = 9.81
    m = 1500
    ro = 1000
    V_volume = 1.4
    A = ro * g * V_volume
    a = 1 - A / (m * g)
    k = 0.5
    
    # Количество шагов
    n_steps = int(T / dt) + 1
    
    # Инициализация массивов
    time = np.linspace(0, T, n_steps)
    
    # Матрица состояния: 8 переменных
    state = np.zeros((n_steps, 8))
    controls = np.zeros((n_steps, 2))  # n_x, n_y
    
    # Начальные условия
    state[0] = [1.0, 0.0, 0.0, 0.0,  # V, theta, x, y
                p0[0], p0[1], p0[2], p0[3]]  # p_V, p_theta, p_x, p_y
    
    # Функция правых частей
    def derivatives(s):
        V, theta, x, y, pV, ptheta, px, py = s
        
        # Защита от деления на 0
        V_safe = max(abs(V), 0.1)
        
        # 1. Вычисляем оптимальное управление
        nx_ideal = k**2 * pV * a * g
        ny_ideal = k**2 * ptheta * a * g / V_safe
        
        # # 2. Применяем ограничения
        # nx = np.clip(nx_ideal, -1.0, 1.0)
        # ny = np.clip(ny_ideal, -3.0, 3.0)
        
        # 3. Вычисляем производные фазовых переменных
        dV = a * g * (nx_ideal - np.sin(theta))
        dtheta = a * g / V_safe * (ny_ideal - np.cos(theta))
        dx = V * np.cos(theta)
        dy = V * np.sin(theta)
        
        # 4. Вычисляем производные сопряженных переменных
        dpx = 0
        dpy = 0
        
        # p_V производная
        dH_dV = px * np.cos(theta) + py * np.sin(theta) - \
                ptheta * (a * g / (V_safe**2)) * (ny_ideal - np.cos(theta))
        dpV = -dH_dV
        
        # p_theta производная
        dH_dtheta = pV * a * g * (-np.cos(theta)) + \
                    ptheta * (a * g / V_safe) * np.sin(theta) + \
                    px * V * (-np.sin(theta)) + py * V * np.cos(theta)
        dptheta = -dH_dtheta
        
        return np.array([dV, dtheta, dx, dy, dpV, dptheta, dpx, dpy]), (nx_ideal, ny_ideal)
    
    # Интегрирование методом Эйлера (можно заменить на РК4)
    for i in range(n_steps - 1):
        derivs, control = derivatives(state[i])
        state[i + 1] = state[i] + derivs * dt
        controls[i] = control
    
    # Последнее управление
    _, last_control = derivatives(state[-1])
    controls[-1] = last_control
    
    return time, state, controls

def simplex_optimization(x_target, y_target, T):
    # Целевая функция для минимизации
    def cost_function(p0):
        # Интегрируем систему с текущими p0
        time, state, _ = integrate_optimal_system(p0, T)
        
        # Конечное состояние
        final = state[-1]
        
        # Ошибки
        pos_error = (final[2] - x_target)**2 + (final[3] - y_target)**2
        adj_error = final[4]**2 + final[5]**2  # p_V(T)^2 + p_theta(T)^2
        
        return pos_error + adj_error
    
    # Начальное предположение
    p0_guess = np.array([0.1, 0.1, 0.0, 0.0])
    
    # Оптимизация симплекс-методом
    result = minimize(cost_function, p0_guess, 
                     method='Nelder-Mead',
                     options={'maxiter': 1500, 
                              'xatol': 1e-4,
                              'fatol': 1e-4,
                              'disp': True})
    
    return result.x, result.success

def plot_optimal_results_simple(time, state, controls, x_target, y_target):
    
    # 1. Скорость V(t)
    print_graphics(time, state[:, 0], 
                   'Оптимальное управление: Скорость V(t)', 
                   'Время, с', 'Скорость, м/с')
    
    # 2. Угол тангажа θ(t)
    print_graphics(time, state[:, 1], 
                   'Оптимальное управление: Угол тангажа θ(t)', 
                   'Время, с', 'Угол, рад')
    
    # 3. Координата x(t)
    print_graphics(time, state[:, 2], 
                   'Оптимальное управление: Координата x(t)', 
                   'Время, с', 'x, м')
    
    # 4. Координата y(t)
    print_graphics(time, state[:, 3], 
                   'Оптимальное управление: Координата y(t)', 
                   'Время, с', 'y, м')
    
    # 5. Управление n_x(t)
    print_graphics(time, controls[:, 0], 
                   'Оптимальное управление: n_x(t)', 
                   'Время, с', 'n_x')
    
    # 6. Управление n_y(t)
    print_graphics(time, controls[:, 1], 
                   'Оптимальное управление: n_y(t)', 
                   'Время, с', 'n_y')
    
    # 7. Траектория движения y(x)
    plt.figure(figsize=(10, 6))
    plt.plot(state[:, 2], state[:, 3], 'b-', linewidth=2)
    plt.plot(x_target, y_target, 'ro', markersize=10, label='Цель')
    plt.plot(state[0, 2], state[0, 3], 'go', markersize=10, label='Старт')
    plt.title('Оптимальное управление: Траектория движения')
    plt.xlabel('x, м')
    plt.ylabel('y, м')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.tight_layout()
    plt.show()
    
    # # 8. Сопряженные переменные p_V(t) и p_theta(t)
    # plt.figure(figsize=(10, 6))
    # plt.plot(time, state[:, 4], 'g-', linewidth=2, label='p_V')
    # plt.plot(time, state[:, 5], 'm-', linewidth=2, label='p_θ')
    # plt.title('Оптимальное управление: Сопряженные переменные')
    # plt.xlabel('Время, с')
    # plt.ylabel('Значение')
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()
    
    # # 9. Сравнение n_x и n_y на одном графике
    # plt.figure(figsize=(10, 6))
    # plt.plot(time, controls[:, 0], 'r-', linewidth=2, label='n_x')
    # plt.plot(time, controls[:, 1], 'b-', linewidth=2, label='n_y')
    # plt.title('Оптимальное управление: Сравнение n_x и n_y')
    # plt.xlabel('Время, с')
    # plt.ylabel('Значение')
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()

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

        fixed_control()
        
        # ====================ШАГ 2====================
        
        optimal_control()
        
