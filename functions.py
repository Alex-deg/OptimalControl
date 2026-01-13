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