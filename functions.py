import os

def additional_test(array, message):
    for func in array:
        if not func(message):
            return False
    return True

def get_valid_input(type, message, arr_of_checking_functions = []):
    while True:
        try:
            if type == 'float' and additional_test(arr_of_checking_functions, message):
                return float(message)
            elif type == 'int' and additional_test(arr_of_checking_functions, message):
                return int(message)
            elif type == 'char' and additional_test(arr_of_checking_functions, message):
                return message
        except ValueError as err:
            if os.name == 'nt':
                os.system('cls')
            elif os.name == 'posix':
                os.system('clear')
            if len(str(err)) > 0:
                print(str(err))
            else:
                print(f"Неверный ввод! Введенное выражение не соответствует типу {type} или не прошло проверку")
            print('Введите значение заново:')
            message = input()


# def print_steps(x_values, y_values, title, x_label, y_label):
#     plt.figure(figsize=(10, 6))
#     plt.step(x_values, y_values, where='post')
#     plt.title(title)
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)
#     plt.tight_layout()
#     plt.show()

# def print_graphics(x_values, y_values, title, x_label, y_label):

#     plt.figure(figsize=(10, 6))
#     plt.plot(x_values, y_values, 'o-', linewidth=2, markersize=8)
#     plt.title(title)
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)
#     plt.tight_layout()
#     plt.show()

# def print_histogram(x_values, y_values, title, x_label, y_label):

#     plt.figure(figsize=(12, 7))

#     bin_width = 0.02

#     plt.bar(x_values, y_values, width=bin_width, 
#             alpha=0.7, color='skyblue',
#             align='center', linewidth=0.5)

#     x_min = (min(x_values) - bin_width/2) * 1.1
#     x_max = (max(x_values) + bin_width/2) * 1.1
#     y_max = max(y_values) * 1.1  
    
#     plt.xlim(x_min, x_max)
#     plt.ylim(0, y_max)
    
#     plt.title(title)
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)
#     plt.grid(True, alpha=0.3, axis='y')
#     plt.tight_layout()
#     plt.show()