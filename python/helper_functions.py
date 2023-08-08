import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from mpmath import sqrt as mpsqrt, mpf, mpc

def commutation(matrix_1, matrix_2):
    '''This function returns the commutation of two matrices given as numpy matrices'''
    return matrix_1@matrix_2 - matrix_2@matrix_1

def give_random_matrix(dim):
    '''This function gives a random Hermitian matrix of dimension dim. So that the first element i.e. element on the first row and first column is zero.'''
    matrix = zero_matrix((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            if i==j and i!= 0:
                matrix[i,j] = (np.random.normal())
            elif i!=j:
                #This was mpc
                random_complex = complex(np.random.normal(), np.random.normal())
                #random_complex = mpc(complex(np.random.normal(), np.random.normal()))
                matrix[i,j] = random_complex
                matrix[j,i] = random_complex.conjugate()
    return matrix

def real_roots(polynomial):
    '''This function returns the roots of a polynomial that are real.'''
    roots = np.roots(polynomial)
    real_roots = []
    for root in roots:
        if np.imag(root) == 0:
            real_roots.append(np.real(root))
    return real_roots

def make_matrix(matrix, size):
    '''This function returns a matrix of given size, with the elements of the matrix given as an argument. The matrix argument must be a list of lists.'''
    return np.matrix(matrix.reshape(size[0], size[1]))

def zero_matrix(size):
    '''This function returns a zero matrix of given size. The size argument must be a tuple of the form (n, m) where n is the number of rows and m is the number of columns.'''
    #This was mpc
    return np.matrix([[complex(0,0) for i in range(size[1])] for j in range(size[0])])
    #return np.matrix([[mpc(complex(0,0)) for i in range(size[1])] for j in range(size[0])])

def kronecker(i, j):
    '''Kronecker delta function returns 1 if the two arguments are equal, return 0 if the two arguments are not equal to each other'''
    return 1 if i==j else 0

def mpc_norm(vector):
    '''This function returns the norm of a vector given as a numpy matrix'''
    norm = (vector@(vector.conjugate().T))[0,0]
    if norm.imag > 0.001:
        return mpsqrt(norm)
    else:
        return np.sqrt(float(norm.real))

def SU2_generators(N):
    '''This function gives the SU(2) generators for arbitriary spin. Spin value is taken as the argument, and gives the three generators as a list of matrices'''
    L_x = []
    L_y = []
    L_z = []
    for k in range(1, int(2*N+1)+1):
        row_x = []
        row_y = []
        row_z = []
        for l in range(1, int(2*N+1)+1):
            row_x.append(np.sqrt((2*N+1-k)*k)*kronecker(k+1,l)/2 + np.sqrt((2*N+1-l)*l)*kronecker(l+1,k)/2)
            row_y.append(-(np.sqrt((2*N+1-k)*k)*kronecker(k+1,l))*1j/2 + (np.sqrt((2*N+1-l)*l)*kronecker(l+1,k))*1j/2)
            row_z.append(((N+1-k)*kronecker(k,l)))
        L_x.append(row_x)
        L_y.append(row_y)
        L_z.append(row_z)
    return [np.array(L_x), np.array(L_y), np.array(L_z)]

def jacobian(func, point, delta=1e-3):
    '''This function returns the Jacobian matrix (as a 2D array given as a numpy matrix) of the function func at the point point. 
    The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)]. both func and point must be numpy matrices.'''
    n = np.size(point)
    result = zero_matrix((n, n))

    for i in range(n):
        for j in range(n):
            result[j, i] = (func(point + np.matrix([0.0 if k!=j else delta for k in range(n)]))[0, i] - func(point - np.matrix([0.0 if k!=j else delta for k in range(n)]))[0, i])/(2*delta)
    return np.matrix(result)

def combined_system(func, point):
    '''This function returns the combined system of the function func and the Jacobian of the function func at the point point. 
    The function func must be a vector function of the form f(x) = [f_1(x), f_2(x), ..., f_n(x)]. both func and point must be numpy matrices.
    matrix should be returned as a flattened np matrix.'''
    size = np.size(func(point)) #dimension of the original system.
    phi = point[0, size:] #the phi matrix which is the last elements of the point matrix.
    phi = make_matrix(phi, (size, size)) #reshape the phi matrix to a square matrix.
    jacobian_matrix = jacobian(func, point[0, :size]) #the jacobian matrix of the original system.
    phi = (phi@jacobian_matrix).flatten() #the new phi matrix, that is evolved with the jacobian matrix, and flattened.
    result = np.concatenate((func(point), phi), axis=1) #return the combined system.
    return result

def split_real_complex(x, dim):
    '''This function splits the real and the complex parts of a list of complex numbers and returns as one list of real numbers.'''
    result = []
    iter_no = int(len(x)/(dim**2))
    for i in range(iter_no):
        reals = []
        imags = []
        for j in range(dim*dim):
            reals.append(np.real(x[i*dim+j]))
            imags.append(np.imag(x[i*dim+j]))
        result.extend(reals)
        result.extend(imags)
    result.extend(np.identity(len(result)).flatten().tolist())
    return result