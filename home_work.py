from math import e, cos, sin
import time
import matplotlib.pyplot as plt
from functions import get_valid_input, additional_test
import numpy as np
from scipy.optimize import minimize, root
from scipy.integrate import solve_ivp

g = 9.81
m = 1500
ro = 1000
V_volume = 1.4
A = ro * g * V_volume
a = 1 - A / (m * g)
k = 0.5

# функция вывода графиков +
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

# first step +
def fixed_control():
    # вводим вектор управления - веткор нагрузок [nx, ny]
    print(f'Введите значения параметров вектора управления для эксперимента через запятую')
    print('nx, ny соответственно')
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
    # заголовки для графиков V, theta, x, y соответственно
    titles_for_graphics = [['График изменения скорости V', 'T, с', 'м/с'],
                           ['График изменения угла тангажа', 'T, c', 'радианы'],
                           ['График изменения координаты по оси OX', 'T, c', ''],
                           ['График изменения координаты по оси OY', 'T, c', '']]
    # цикл вывода графиков интересующих нас значений - V, theta, x, y
    x_axis = np.arange(0, T, 0.01)
    for i in range(len(list_of_values)):
        main_title, x_title, y_title = titles_for_graphics[i]
        print_graphics(x_axis, list_of_values[i], main_title, x_title, y_title)
    print_graphics(list_of_values[2], list_of_values[3], 'График зависимости y от x', 'x', 'y')

def method_processing(method_name, method_func, x_f, y_f, T, p0_guess=[]):

    result = {}

    start_time = time.time()
    if method_name == 'Newtone Hybrid':
        p0_optimal, success, iterations = method_func(x_f, y_f, T, p0_guess) 
    else:
        p0_optimal, success, iterations = method_func(x_f, y_f, T)
    end_time = time.time()
    
    if success:
        _, state, _ = integrate_optimal_system(p0_optimal, T, method='RK45')
        final = state[-1]
        
        pos_error = np.sqrt((final[2] - x_f) ** 2 + (final[3] - y_f) ** 2)
        adj_error = np.sqrt(final[4] ** 2 + final[5] ** 2)
        
        print(f"  Время вычисления: {end_time - start_time:.3f} с")
        print(f"  Итераций: {iterations}")
        print(f"  Ошибка позиции: {pos_error:.6e}")
        print(f"  Ошибка сопряженных: {adj_error:.6e}")
        print(f"  p0_optimal: {p0_optimal}")
        print(f"  Конечная позиция: x = {final[2]}; y = {final[3]}")
        
        result = {
            'p0': p0_optimal,
            'time': end_time - start_time,
            'iterations': iterations,
            'pos_error': pos_error,
            'adj_error': adj_error
        }

        # Вывод графиков
        print('\nВывести графики? (y/n): ', end='')
        if input().lower() == 'y':
            time_vals, state_vals, controls_vals = integrate_optimal_system(p0_optimal, T, method='RK45')
            plot_optimal_results_simple(time_vals, state_vals, controls_vals, x_f, y_f)

    else:
        print(f"  Метод {method_name} не сошелся!")
        result = None

    return result

# second step
def optimal_control():
    print(f'Введите координаты точки, в которую хотите попасть:')
    print('x_f, y_f = ', end='')
    x_f, y_f = map(float, input().split(', '))
    print(f'Введите время за которое хотите добраться из [0; 0] в [{x_f}; {y_f}]')
    T = float(input())
    
    print("\n" + "="*40)
    print("1. СРАВНЕНИЕ МЕТОДОВ")
    print("="*40)
    
    # Параметры для сравнения
    methods = [
        ('СИМПЛЕКС (Nelder-Mead)', simplex_optimization),
        ('НЬЮТОН', newton_method_optimization)
    ]
    
    results = {}
    
    for method_name, method_func in methods:
        print(f"\n{method_name}:")
        print("-" * 30)
        results[method_name] = method_processing(method_name, method_func, x_f, y_f, T)
    
    # 5. Гибридный подход: Симплекс → Ньютон
    print("\n" + "="*40)
    print("ГИБРИДНЫЙ ПОДХОД: Симплекс -> Ньютон")
    print("="*40)
    
    if 'СИМПЛЕКС (Nelder-Mead)' in results and results['СИМПЛЕКС (Nelder-Mead)'] is not None:
        simplex_result = results['СИМПЛЕКС (Nelder-Mead)']
        p0_simplex = simplex_result['p0']
        results['Newtone Hybrid'] = method_processing('Newtone Hybrid', newton_method_optimization, x_f, y_f, T, p0_simplex)

# Функция правых частей
def derivatives(t, s):
    V, theta, x, y, pV, ptheta, px, py = s
    
    # Защита от деления на 0
    V_safe = max(abs(V), 0.1)
    
    # 1. Вычисляем оптимальное управление
    nx = k**2 * pV * a * g
    ny = k**2 * ptheta * a * g / V_safe
    
    # # 2. Применяем ограничения
    # nx = np.clip(nx_ideal, -1.0, 1.0)
    # ny = np.clip(ny_ideal, -3.0, 3.0)
    
    # 3. Вычисляем производные фазовых переменных
    dV = a * g * (nx - np.sin(theta))
    dtheta = a * g / V_safe * (ny - np.cos(theta))
    dx = V * np.cos(theta)
    dy = V * np.sin(theta)
    
    # 4. Вычисляем производные сопряженных переменных
    dpx = 0
    dpy = 0
    
    # p_V производная
    dH_dV = px * np.cos(theta) + py * np.sin(theta) - \
            ptheta * (a * g / (V_safe**2)) * (ny - np.cos(theta))
    dpV = -dH_dV
    
    # p_theta производная
    dH_dtheta = pV * a * g * (-np.cos(theta)) + \
                ptheta * (a * g / V_safe) * np.sin(theta) + \
                px * V * (-np.sin(theta)) + py * V * np.cos(theta)
    dptheta = -dH_dtheta
    
    return np.array([dV, dtheta, dx, dy, dpV, dptheta, dpx, dpy]), (nx, ny)

def integrate_optimal_system(p0, T, method='RK45', dt=0.01, rtol=1e-6, atol=1e-9):
    
    def system_dynamics(t, z):
        V, theta, x, y, pV, ptheta, px, py = z
        
        # Защита от деления на ноль
        V_safe = max(abs(V), 1e-6)
        
        # 1. Вычисление оптимального управления по принципу максимума
        nx = k**2 * pV * a * g
        ny = k**2 * ptheta * a * g / V_safe
        
        # 2. Производные фазовых переменных
        dV = a * g * (nx - np.sin(theta))
        dtheta = a * g / V_safe * (ny - np.cos(theta))
        dx = V * np.cos(theta)
        dy = V * np.sin(theta)
        
        # 3. Производные сопряженных переменных
        dpx = 0
        dpy = 0
        
        # Производная для pV
        dH_dV = px * np.cos(theta) + py * np.sin(theta) - \
                ptheta * (a * g / (V_safe**2)) * (ny - np.cos(theta))
        dpV = -dH_dV
        
        # Производная для ptheta
        dH_dtheta = pV * a * g * (-np.cos(theta)) + \
                   ptheta * (a * g / V_safe) * np.sin(theta) + \
                   px * V * (-np.sin(theta)) + py * V * np.cos(theta)
        dptheta = -dH_dtheta
        
        return [dV, dtheta, dx, dy, dpV, dptheta, dpx, dpy]
    
    # Начальные условия
    V0 = 1.0
    theta0 = 0.0
    x0 = 0.0
    y0 = 0.0
    z0 = [V0, theta0, x0, y0, p0[0], p0[1], p0[2], p0[3]]
    
    # Интегрирование системы
    t_eval = np.arange(0, T + dt, dt)
    sol = solve_ivp(system_dynamics, [0, T], z0, 
                    method=method,
                    t_eval=t_eval,
                    rtol=rtol, 
                    atol=atol)
    
    # Извлекаем результаты
    time = sol.t
    state = sol.y.T  # транспонируем для удобства: (n_steps, 8)
    
    # Вычисляем управления для каждой точки
    n_steps = len(time)
    controls = np.zeros((n_steps, 2))
    for i in range(n_steps):
        V, theta, _, _, pV, ptheta, _, _ = state[i]
        V_safe = max(abs(V), 1e-6)
        controls[i, 0] = k**2 * pV * a * g
        controls[i, 1] = k**2 * ptheta * a * g / V_safe
    
    return time, state, controls

def plot_optimal_results_simple(time, state, controls, x_target, y_target):
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 12))
    fig.suptitle('Результаты оптимального управления', fontsize=16)
    
    # 1. Скорость V(t)
    axes[0, 0].plot(time, state[:, 0], 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Время, с')
    axes[0, 0].set_ylabel('Скорость, м/с')
    axes[0, 0].set_title('Скорость V(t)')
    axes[0, 0].grid(True)
    
    # 2. Угол тангажа θ(t)
    axes[0, 1].plot(time, np.degrees(state[:, 1]), 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Время, с')
    axes[0, 1].set_ylabel('Угол, град')
    axes[0, 1].set_title('Угол тангажа θ(t)')
    axes[0, 1].grid(True)
    
    # 3. Координата x(t)
    axes[0, 2].plot(time, state[:, 2], 'g-', linewidth=2)
    axes[0, 2].plot([time[-1]], [x_target], 'ro', markersize=10)
    axes[0, 2].set_xlabel('Время, с')
    axes[0, 2].set_ylabel('x, м')
    axes[0, 2].set_title(f'Координата x(t) (цель: {x_target} м)')
    axes[0, 2].grid(True)
    
    # 4. Координата y(t)
    axes[1, 0].plot(time, state[:, 3], 'm-', linewidth=2)
    axes[1, 0].plot([time[-1]], [y_target], 'ro', markersize=10)
    axes[1, 0].set_xlabel('Время, с')
    axes[1, 0].set_ylabel('y, м')
    axes[1, 0].set_title(f'Координата y(t) (цель: {y_target} м)')
    axes[1, 0].grid(True)
    
    # 5. Управление n_x(t)
    axes[1, 1].plot(time, controls[:, 0], 'c-', linewidth=2)
    axes[1, 1].set_xlabel('Время, с')
    axes[1, 1].set_ylabel('n_x')
    axes[1, 1].set_title('Управление n_x(t)')
    axes[1, 1].grid(True)
    
    # 6. Управление n_y(t)
    axes[1, 2].plot(time, controls[:, 1], 'y-', linewidth=2)
    axes[1, 2].set_xlabel('Время, с')
    axes[1, 2].set_ylabel('n_y')
    axes[1, 2].set_title('Управление n_y(t)')
    axes[1, 2].grid(True)
    
    # 7. Траектория движения
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
    
def newton_method_optimization(x_target, y_target, T, p0_guess=None, method='hybr'):
    
    if p0_guess is None:
        p0_guess = np.array([0.1, 0.1, 0.0, 0.0])
    
    def boundary_conditions(p0):

        _, state, _ = integrate_optimal_system(p0, T, method='RK45', rtol=1e-6, atol=1e-9)
        
        # Конечное состояние
        zT = state[-1]
        
        # Невязки граничных условий
        F1 = zT[2] - x_target  # x(T) - x_target
        F2 = zT[3] - y_target  # y(T) - y_target
        F3 = zT[4]             # pV(T) (должно быть 0)
        F4 = zT[5]             # ptheta(T) (должно быть 0)
        
        return np.array([F1, F2, F3, F4])
    
    result = root(boundary_conditions, p0_guess, 
                  method=method, 
                  tol=1e-8,     
                  options={
                      'maxfev': 2000, 
                      'xtol': 1e-8, 
                      'eps': 1e-8,
                      'factor': 0.1
                  })
    
    return result.x, result.success, result.nfev

def simplex_optimization(x_target, y_target, T, method='Nelder-Mead'):
    
    def cost_function(p0):
        _, state, _ = integrate_optimal_system(p0, T, method='RK45', rtol=1e-6, atol=1e-9)
        
        # Конечное состояние
        final = state[-1]
        
        # Ошибки с весами
        pos_error = (final[2] - x_target)**2 + (final[3] - y_target)**2
        adj_error = final[4]**2 + final[5]**2  # p_V(T)^2 + p_theta(T)^2
        
        return pos_error + 0.1 * adj_error  # весовой коэффициент
    
    # Начальное приближение
    p0_guess = np.array([0.1, 0.1, 0.0, 0.0])
    
    # Оптимизация с теми же параметрами точности
    result = minimize(cost_function, p0_guess, 
                     method=method,
                     options={
                         'maxiter': 2000, 
                         'xatol': 1e-8,   
                         'fatol': 1e-8,    
                         'disp': False
                     })
    
    return result.x, result.success, result.nit

if __name__ == "__main__":
    
    # считываем общее количество экспериментов, в которых будем проверять заисимость V, theta, x и y от [nx, ny, nz]
    print('Введите количество экспериментов\nexperiments = ', end='')  
    experiments = input()
    experiments = get_valid_input('int', experiments, [])

    # основной цикл
    for i in range(experiments):
        # ====================ШАГ 1====================
        #fixed_control()     
        # ====================ШАГ 2====================
        optimal_control()