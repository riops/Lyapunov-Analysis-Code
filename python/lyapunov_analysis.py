import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from mpmath import sqrt as mpsqrt, mpf, mpc
from helper_functions import *
from systems import *

def RKStep(func, x, Dt):
    '''This function returns the next point in the Runge-Kutta method. The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)].
    x and t must be numpy matrices. Dt is the step size. It does it for the combined system.'''
    k1 = combined_system(func, x) * Dt
    k2 = combined_system(func, x + k1 * 0.5) * Dt
    k3 = combined_system(func, x + k2 * 0.5) * Dt
    k4 = combined_system(func, x + k3) * Dt
    return x + (k1 + 2*k2 + 2*k3 + k4)/6

def RKiterate(func, x, Dt, number_of_steps, give_all_points=False, print_progress=False):
    '''This function returns the next point in the Runge-Kutta method. The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)].
    x and t must be numpy matrices. Dt is the step size. It does it for the combined system.'''
    time_list = []
    all_points = np.matrix(zero_matrix((number_of_steps + 1, np.size(x))))
    all_points[0, :] = x
    for step in range(number_of_steps):
        print('Current iteration for the integral is {}...'.format(step+1), end='\r') if print_progress else None
        time_list.append(step*Dt)
        x = RKStep(func, x, Dt)
        #print(energy_check(make_matrix(x[0,:4], (2,2))))
        all_points[step + 1, :] = x
    if give_all_points:
        return all_points, time_list
    else:
        return x

def project(vector_1, vector_2):
    '''This function projects vector_1 onto vector_2. Both vectors must be numpy matrices.'''
    return ((vector_1.conjugate()@vector_2.T)/(vector_1@vector_1.conjugate().T))*vector_1

def gram_schmidt(vectors, normalize=False):
    '''This function takes a list of vectors as a np matrix and returns a list of orthogonal vectors as a np matrix. normalize condition normalizes the vectors.
    Function takes the rows as vectors.'''
    number_of_vectors = np.size(vectors, 0)
    orthogonal_vectors = np.matrix(zero_matrix((number_of_vectors, np.size(vectors, 1))))
    if normalize:
        for i in range(number_of_vectors):
            orthogonal_vectors[i, :] = vectors[i, :]
            if mpc_norm(orthogonal_vectors[i, :]) != 0:
                for j in range(i):
                    orthogonal_vectors[i, :] = orthogonal_vectors[i, :] - project(orthogonal_vectors[j, :], vectors[i, :])
            else:
                orthogonal_vectors[i, :] = np.matrix(zero_matrix((1, np.size(vectors, 1))))
        for j in range(number_of_vectors):
            orthogonal_vectors[j, :] = orthogonal_vectors[j, :]/mpc_norm(orthogonal_vectors[j, :])
        return orthogonal_vectors
    else:
        for i in range(number_of_vectors):
            orthogonal_vectors[i, :] = vectors[i, :]
            if mpc_norm(orthogonal_vectors[i, :]) != 0:
                for j in range(i):
                    orthogonal_vectors[i, :] = orthogonal_vectors[i, :] - project(orthogonal_vectors[j, :], vectors[i, :])
            else:
                orthogonal_vectors[i, :] = np.matrix(zero_matrix((1, np.size(vectors, 1))))
        return orthogonal_vectors

def make_float(x):
    '''This function makes a matrix of floats from a matrix of complex numbers.'''
    x = x.tolist() if type(x) == np.matrix else x
    real_matrix = []
    for row in x:
        real_row = []
        for element in row:
            real_row.append((float(element.real)))
        real_matrix.append(real_row)
    return np.matrix(real_matrix)

def largest_lyapunov(func, x, Dt, number_of_steps):
    '''This function gives only the largest lypunov exponent. The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)].'''
    size = np.size(func(x))
    final_point = RKiterate(func, x, Dt, number_of_steps, give_all_points=False)
    phi_final = make_matrix(final_point[0, size:], (size, size))
    random_vector = np.matrix(np.random.rand(size)).T
    lyapunov_vector = phi_final@random_vector
    return np.log(mpc_norm(lyapunov_vector))/(number_of_steps*Dt)

def lyapunov_spectrum(func, x, Dt, number_of_steps, normalization_number, give_all_points=False, print_progress=True):
    '''This function gives the lypunov spectrum. The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)].'''
    size = np.size(func(x)) #We get the size of the system.
    lyapunovs = np.matrix(zero_matrix((normalization_number, size))) #We initialize the lyapunovs matrix. lyapunovs[i, j] is the jth lyapunov exponent of the ith time step.
    logs = np.matrix(zero_matrix((normalization_number, size))) #We initialize the logs matrix. logs[i, j] logarithm of the norm of the jth vector exponent of the ith time step.
    for i in range(normalization_number):
        x = RKiterate(func, x, Dt, number_of_steps) #We integrate the system from the initial point to the final point, with number_of_steps steps, and Dt step size.
        phi_intermediate = make_matrix(x[0, size:], (size, size)) #We convert the phi part of the point to a matrix, so that we can gram-schmidt orthogonalize it.
        phi_intermediate = gram_schmidt(phi_intermediate) #We orthogonalize the phi part of the point.
        print('Current normalization number is {}...'.format(i+1), end='\r') if print_progress else None
        for j in range(size):
            norm = mpc_norm(phi_intermediate[j, :]) #For each vector we calculate the norm.
            logs[i, j] = np.log(norm) #We calculate the log of the norm of the vector.
            lyapunovs[i, j] = (np.sum(logs[:i, j]))/((i+1)*number_of_steps*Dt) #We calculate the lyapunov exponent. We do this by summing the logs of the norms of the vectors, and dividing by the number of steps, how many times we have normalized and the step size.
            phi_intermediate[j,:] = phi_intermediate[j,:]/norm #We normalize the vectors. One by one.
        x[0, size:] = phi_intermediate.flatten() #We convert the matrix back to a vector, and put it in the point, so that it is taken as the initial point for the next iteration.
    if give_all_points:
        return make_float(lyapunovs).T.tolist() #We return the whole lyapunovs matrix, which is the lyapunov spectrum.
    else:
        return lyapunovs[-1, :] #We return only the last row of the lyapunovs matrix, which is the largest lyapunov exponent.

def lyapunov_analysis(energy_function, function, energy_range, energy_check, Dt=0.01, number_of_steps=1, normalization_number=3_000, number_of_points=10, print_progress=True, only_largest=True):
    '''This function calculates the mean largest lyapunov exponent for a given energy range. It also calculates the error of the mean. 
    The function energy_function must be a function that takes the energy as an argument and returns the initial conditions for the system.
    The function function must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)]. Which is the function we integrate over.
    The energy_range must be a list of the form [min_energy, max_energy, increment]. The energy range is the range of energies we want to calculate the lyapunov spectrum for.
    The energy_check is the energy we want to check that the initial conditions are correct for.
    The Dt is the step size of the integration.
    The number_of_steps is the number of steps we do not normalize the vectors for.
    The normalization_number is the number of times we normalize the vectors.
    The number_of_points per energy is the number of points we want to calculate the lyapunov spectrum for.
    The get_energy_from_file is a boolean that tells us if we want to get the energies from a file or not.
    The print_progress is a boolean that tells us if we want to print the progress of the calculation or not.
    The only_largest is a boolean that tells us if we want to only calculate the largest lyapunov exponent or not.'''
    all_lyapunovs = [] #We initialize the list of all lyapunovs.
    max_lyapunovs = [] #We initialize the list of the largest lyapunov exponents.
    average_max_lyapunovs = [] #We initialize the list of the average of the largest lyapunov exponents.
    errors = [] #We initialize the list of the errors of the largest lyapunov exponents.
    for energy in energy_range: #We iterate over the energy range.
        initial_coditions = [] #We initialize the list of initial conditions.
        point = 0 #We initialize the number of points.
        while point < number_of_points: #We iterate until we have enough points.
            initial_condition = energy_function(energy) #We get an initial condition with the given energy.
            energy_difference = abs(energy_check(initial_condition) - energy) #We calculate the difference between the energy of the initial condition and the given energy.
            if energy_difference < 0.01: #We check if the energy of the initial condition is close enough to the given energy.
                print('Initial condition for energy {} and point {} is correct.'.format(energy, point+1))
                initial_coditions.append(initial_condition) #If it is, we add it to the list of initial conditions.
                point += 1 #We increase the number of points.
        inputs = [(function, initial_coditions[i], Dt, number_of_steps, normalization_number) for i in range(number_of_points)] #We create the inputs for the lyapunov_spectrum function.
        with Pool() as p: #We use multiprocessing to calculate the lyapunov spectrum of each initial condition.
            lyapunovs = p.starmap(lyapunov_spectrum, inputs) #We calculate the lyapunov spectrum of each initial condition.
            all_lyapunovs.append(lyapunovs) #We add the lyapunov spectrum of the initial conditions to the list of all lyapunovs.
            max_lyapunovs_for_energy = []
            for lyapunov in lyapunovs:
                max_lyapunov = max(make_float(lyapunov).tolist()[0]) #We calculate the largest lyapunov exponent.
                max_lyapunovs.append(max_lyapunov) #We add the largest lyapunov exponent to the list of the largest lyapunov exponents.
                max_lyapunovs_for_energy.append(max_lyapunov) #We add the largest lyapunov exponent to the list of the largest lyapunov exponents for the given energy and given points.
        average_max_lyapunov = sum(max_lyapunovs_for_energy)/len(max_lyapunovs_for_energy) #We calculate the average of the largest lyapunov exponents.
        average_max_lyapunovs.append(average_max_lyapunov) #We add the average of the largest lyapunov exponents to the list of the average of the largest lyapunov exponents.
        standart_deviation = np.std(max_lyapunovs_for_energy) #We calculate the standart deviation of the largest lyapunov exponents.
        error = standart_deviation/np.sqrt(len(max_lyapunovs_for_energy)) #We calculate the error of the largest lyapunov exponents.
        errors.append(error) #We add the error of the largest lyapunov exponents to the list of the errors of the largest lyapunov exponents.
        print('For the energy value {} the average largest lyapunov exponent is {} with error {}'.format(energy, average_max_lyapunov, error)) if print_progress else None #We print the average of the largest lyapunov exponents and the error of the largest lyapunov exponents.
    if only_largest: #If we only want to return the largest lyapunov exponents, we do this.
        return average_max_lyapunovs, errors #We return the average of the largest lyapunov exponents and the error of the largest lyapunov exponents.
    return all_lyapunovs #We return the list of all lyapunovs.

def largest_lyapunov_analysis(energy_function, function, energy_range, Dt=0.01, number_of_steps = 1_000, number_of_points=1, print_progress=True):
    '''This function calculates the largest lyapunov exponent for a given energy range.
    The energy_function is the function that gives the initial conditions for a given energy.
    The function is the function that we want to calculate the largest lyapunov exponent for.
    The energy_range is the range of energies we want to calculate the largest lyapunov exponent for.
    The Dt is the time step we want to use.
    The number_of_steps is the number of steps we want to use.
    The number_of_points is the number of points per energy we want to use.
    The print_progress is a boolean that tells us if we want to print the progress or not.'''
    all_largest_lyaps = []
    for energy in energy_range:
        for i in range(number_of_points):
            print('Starting calculation for energy {} and point {}'.format(energy, i+1)) if print_progress else None
            largest_lyapunovs = largest_lyapunov(function, energy_function(energy), Dt, number_of_steps)
            print('Largest lyapunov exponent for energy {} calculated. It is {}'.format(energy, largest_lyapunovs)) if print_progress else None
            all_largest_lyaps.append(largest_lyapunovs)
    return all_largest_lyaps

def plot_lyapunov_spectrum(lyapunov_spectrum, time_list, title=None, x_label=None, y_label=None, save=False, save_name=None, show=True):
    '''This function plots the lyapunov spectrum.
    The lyapunov_spectrum is the lyapunov spectrum we want to plot.
    The title is the title of the plot.
    The x_label is the label of the x axis.
    The y_label is the label of the y axis.
    The save is a boolean that tells us if we want to save the plot or not.
    The save_name is the name of the file we want to save the plot as.
    The show is a boolean that tells us if we want to show the plot or not.'''
    for i in range(len(lyapunov_spectrum)):
        plt.plot(time_list, lyapunov_spectrum[i])
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if save:
        plt.savefig(save_name)
    elif show:
        return plt    