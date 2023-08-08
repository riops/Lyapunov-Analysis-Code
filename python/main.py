from lyapunov_analysis import *
import functools as ft
import inspect
import os
import time
from datetime import datetime

def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }

def change_default_args(func, dict_of_args):
    for arg in dict_of_args:
        func = ft.partial(func, **{arg: dict_of_args[arg]})
    return func

def select_system(system_select_input=''):
    system_dictionary = {'toy_model': [toy_model, initial_condition_toy_model, hamiltonian_toy_model], 
                         'new_model': [new_model, initial_condition_new_model, hamiltonian_new_model],
                         'original_model': [F, initial_condition_F, hamiltonian_F],
                         'original_model_exact': [G, initial_condition_G, hamiltonian_G]}
    if system_select_input == '':
        system_select_input = input("Please write the name of the system you want to use: ")
        system_select_input = system_select_input.replace(' ', '_').lower()
    else:
        system_select_input = system_select_input.replace(' ', '_').lower()
    
    if system_select_input == 'help' or system_select_input == 'h':
        print('Names of the available systems are: ')
        for system_name in system_dictionary:
            print(system_name)
        system_select_input = input("Please write the name of the system you want to use: ")
        return select_system(system_select_input)
    elif system_select_input == 'exit':
        return False
    else:
        try:
            system = system_dictionary[system_select_input][0]
            initial_condition_func = system_dictionary[system_select_input][1]
            hamiltonian = system_dictionary[system_select_input][2]
            return system, initial_condition_func, hamiltonian, system_select_input
        except KeyError:
            print("System not found. Please try again.")
            return select_system()

def change_default_arguments(system_functions):
    default_argument_dict = get_default_args(system_functions[0])
    print('System you have selected is: {}, and it has the following deafult arguments {}'.format(system_functions[3], default_argument_dict))
    default_argument_change_option = input("Do you want to change the default arguments? (y/n): ").lower()

    selected_system = system_functions[3]
    if default_argument_change_option == 'y':
        for arg in default_argument_dict:
            new_value = input('Please write the new value for the argument {} (Press enter to keep the default value): '.format(arg))
            if new_value != '':
                default_argument_dict[arg] = float(new_value)
            system = change_default_args(system_functions[0], default_argument_dict)
            initial_condition_func = change_default_args(system_functions[1], default_argument_dict)
            hamiltonian = change_default_args(system_functions[2], default_argument_dict)            
        return change_default_arguments((system, initial_condition_func, hamiltonian, selected_system))
    elif default_argument_change_option == 'n':
        return system_functions
    elif default_argument_change_option == 'exit':
        return False
    else:
        print('Please write y, n or exit')
        return change_default_arguments(system)

def change_default_parameters():
    parameter_dictionary = {
                            'iter_no': 2_000,
                            'Dt': 0.01,
                            'numer_of_iters_before_normalizing': 10,
                            'energy_tolerance': 0.1,
                            'lowest_energy': 100, 
                            'highest_energy': 1100, 
                            'energy_step': 100,
                            'energy': 100,
                            'N': 10
                            }
    default_parameter_change_option = input("Do you want to change the default parameters? (y/n): ").lower()
    if default_parameter_change_option == 'help' or default_parameter_change_option == 'h':
        print('Default parameters are: ')
        for parameter in parameter_dictionary:
            print('{}: {}'.format(parameter.replace('_', ' '), parameter_dictionary[parameter]))
        return change_default_parameters()
    elif default_parameter_change_option == 'y':
        for parameter in parameter_dictionary:
            new_value = input('Please write the new value for the parameter {} (Press enter to keep the default value): '.format(parameter.replace('_', ' ')))
            if new_value != '':
                parameter_dictionary[parameter] = float(new_value)
        return parameter_dictionary
    elif default_parameter_change_option == 'n':
        return parameter_dictionary
    elif default_parameter_change_option == 'exit':
        return False
    else:
        print('Please write y, n or exit')
        return change_default_parameters()

def analyze_system(system_functions, parameters):
    iter_no = int(parameters['iter_no'])
    Dt = parameters['Dt']
    numer_of_iters_before_normalizing = int(parameters['numer_of_iters_before_normalizing'])
    energy_tolerance = parameters['energy_tolerance']
    energy_range = range(int(parameters['lowest_energy']), int(parameters['highest_energy']), int(parameters['energy_step']))
    energy = int(parameters['energy'])
    N = int(parameters['N'])

    system = system_functions[0]
    initial_condition_func = system_functions[1]
    hamiltonian = system_functions[2]
    selected_system = system_functions[3]

    analysis_mode = input("Please write the analysis mode (1 for plotting the spectrum, 2 for getting the largest lyapunov analysis): ").lower()
    if analysis_mode == '1':
        start = time.time()
        while True:
            initial_condition = initial_condition_func(energy)
            if abs(hamiltonian(initial_condition) - energy) < energy_tolerance:
                break
        lyapunov_spectrum_list = lyapunov_spectrum(system, initial_condition, Dt, numer_of_iters_before_normalizing, iter_no, give_all_points=True)
        time_list =  [i*Dt*numer_of_iters_before_normalizing for i in range(iter_no)]
        now = datetime.now()
        path = os.path.join(os.getcwd(), 'Lyapunov_spectrums_{}'.format(selected_system))
        if not os.path.exists(path):
            os.mkdir(path)
        dt_string = now.strftime('%d_%m_%Y_%H_%M_%S')
        save_name = os.path.join(path, 'lyapunov_spectrum_dim={}_energy={}_date={}.png'.format(2, energy, dt_string))
        plot_lyapunov_spectrum(lyapunov_spectrum_list, time_list, show=True, save=True, save_name=save_name)
        print('Time taken is: {}'.format(time.time() - start))
        return True
    elif analysis_mode == '2':
        start = time.time()
        lyaps = lyapunov_analysis(initial_condition_func, system, energy_range, hamiltonian, Dt=Dt, number_of_steps=numer_of_iters_before_normalizing, normalization_number=iter_no)
        mean_largest_lyapunovs = lyaps[0]
        print(mean_largest_lyapunovs)
        error_list = lyaps[1]
        print(error_list)
        print('Time taken is: {}'.format(time.time() - start))  
        return True
    elif analysis_mode == 'exit':
        return False
    else:
        print('Please write 1, 2 or exit')
        return analyze_system(system_functions, parameters)

if __name__ == '__main__':
    continue_program = True
    while continue_program:
        system = select_system()
        if not system:
            continue_program = False
            break
        else:
            system = change_default_arguments(system)
            if not system:
                continue_program = False
                break
            else:
                parameters = change_default_parameters()
                result = analyze_system(system, parameters)
                if not result:
                    continue_program = False
                    break

                continue_program = input("Do you want to continue? (y/n): ").lower()
                if continue_program == 'y':
                    continue_program = True
                elif continue_program == 'n':
                    continue_program = False
                elif continue_program == 'exit':
                    continue_program = False
                    break