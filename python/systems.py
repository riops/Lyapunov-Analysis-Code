import numpy as np
from helper_functions import *
import random
import math
   
def F(x, m1=100, m2=1, dim=2):
    '''This is fourth order model with the non-chaotic behavior i.e. our original model.'''
    position = make_matrix(x[0,:dim**2], (dim, dim))
    momentum = make_matrix(x[0,dim**2:2*dim**2], (dim, dim))

    L = SU2_generators((dim-1)/2)
    L_x, L_y, L_z = L[0], L[1], L[2]

    position_dot = momentum.T.flatten()
    momentum_dot = (commutation(L_x, commutation(L_x, position)) + commutation(L_y, commutation(L_y, position)) + commutation(L_z, commutation(L_z, position))) - m1 * position - (2 * m2) * (position@position@position)
    momentum_dot = momentum_dot.T.flatten()
    result = np.concatenate((position_dot, momentum_dot), axis=1)
    return result

def initial_condition_F(E, m1=100, m2=1, dim=2):
    '''This function gives the initial conditions for a given energy E, for the original model'''
    phi = give_random_matrix(dim)
    Pphi = give_random_matrix(dim)
    L = SU2_generators((dim-1)/2)
    polynomial_pos = []
    polynomial_pos.append(m2/2) #Adding the fourth order term
    polynomial_pos.append(0) #Adding the third order term
    second_order = (m1/2 + (dim**2-1)/4) #Adding the second order term
    for i in range(1, dim):
        second_order += (phi[0,i]*phi[i,0])*(2*m2)
    for j in range(3):
        second_order -= (L[j][0,0]**2)
    polynomial_pos.append(np.real(second_order))
    first_order = 0 #Adding the first order terms
    for k in range(1, dim):
        for l in range(1, dim):
            first_order += (phi[0,k]*phi[k,l]*phi[l,0])*(2*m2)
    for i in range(3):
        for j in range(dim):
            for k in range(dim):
                if k != 0 or j != 0:
                    first_order -= -2*L[i][j,0]*L[i][0,k]*phi[k,j]
    polynomial_pos.append(np.real(first_order))
    zeroth_order = np.trace(phi@phi@phi@phi)*(m2/2) #Adding the zeroth order terms
    zeroth_order += np.trace(phi@phi)*(m1/2)
    for i in range(3):
        zeroth_order += np.trace(commutation(L[i], phi)@commutation(L[i], phi))*(-1/2)
    zeroth_order += -E/2
    polynomial_pos.append(np.real(zeroth_order))
    p_list = real_roots(polynomial_pos)
    zeroth_order_momentum = np.trace(Pphi@Pphi)/2-E/2
    polynomial_mom = [1/2, 0, zeroth_order_momentum]
    p_list_mom = real_roots(polynomial_mom)
    if len(p_list) != 0 and len(p_list_mom) != 0:
        r = p_list[random.randint(0, len(p_list)-1)]
        q = p_list_mom[random.randint(0, len(p_list_mom)-1)]
        #this was mpc
        phi[0,0] = complex(r,0)
        Pphi[0,0] = complex(q,0)
        #Pphi = zero_matrix((dim, dim))
        PPhi = Pphi.flatten()
        phi = phi.flatten()
        phi = np.concatenate((phi,PPhi), axis=1)
        phi = np.concatenate((phi, np.matrix(np.identity(2*dim**2)).flatten()), axis=1)
        return phi
    else:
        return initial_condition_F(E)

def hamiltonian_F(initial_condition, m1=100, m2=-1, dim=2):
    '''This function gives the hamiltonian for a given initial condition for the original model.'''
    momentum = make_matrix(initial_condition[0,dim**2:2*dim**2], (dim, dim))
    position = make_matrix(initial_condition[0,:dim**2], (dim, dim))
    L = SU2_generators((dim-1)/2)
    energy = np.trace(position@position@position@position)*(m2/2)
    energy += np.trace(position@position)*(m1/2)
    energy += np.trace(momentum@momentum)/2
    for i in range(3):
        energy -= np.trace(commutation(L[i], position)@commutation(L[i], position))/2
    return float(energy.real)

def G(x, m1=1, m2=60):
    '''This is the exact solution for the original system with matrix dimension set to 2.'''
    return 2*np.matrix([x[0,4], 
                      x[0,5], 
                      x[0,6], 
                      x[0,7], 
                      -(x[0,0]*m1 + 2*m2*x[0,0]**3 + 6*m2*x[0,0]*(x[0,1]**2+x[0,2]**2+x[0,3]**2)), 
                      -(2*x[0,1] + m1*x[0,1] + m2*x[0,1]*(x[0,1]**2+x[0,2]**2+x[0,3]**2)+ m2*x[0,1]*(6*x[0,0]**2+x[0,1]**2+x[0,2]**2+x[0,3]**2)), 
                      -(2*x[0,2] + m1*x[0,2] + m2*x[0,2]*(x[0,1]**2+x[0,2]**2+x[0,3]**2)+ m2*x[0,2]*(6*x[0,0]**2+x[0,1]**2+x[0,2]**2+x[0,3]**2)), 
                      -(2*x[0,3] + m1*x[0,3] + m2*x[0,3]*(x[0,1]**2+x[0,2]**2+x[0,3]**2)+ m2*x[0,3]*(6*x[0,0]**2+x[0,1]**2+x[0,2]**2+x[0,3]**2))])

def initial_condition_G(E, m1=1, m2=60):
    '''This function gives an initial condition for the exact solution of the original model, for a given energy E.'''
    poly = []
    Pphi = np.matrix([[0 for i in range(2)] for j in range(2)]).flatten()
    phi = np.matrix([[0, np.random.normal()], [np.random.normal(), np.random.normal()]])
    poly.append(m2)
    poly.append(0)
    poly.append(m1+6*m2*(phi[0,1]**2+phi[1,0]**2+phi[1,1]**2))
    poly.append(0)
    poly.append((m1+2)*((phi[0,1]**2+phi[1,0]**2+phi[1,1]**2) + m2*(phi[0,1]**2+phi[1,0]**2+phi[1,1]**2))-2*E)
    p_list = real_roots(poly)
    if len(p_list) != 0:
        r = p_list[random.randint(0, len(p_list)-1)]
        phi[0,0] = r
        phi = phi.flatten()
        phi = np.concatenate((phi,Pphi), axis=1)
        phi = np.concatenate((phi, np.matrix(np.identity(8)).flatten()), axis=1)
        return phi.tolist()[0]
    else:
        return initial_condition_G(E)

def hamiltonian_G(x, m1=1, m2=60):
    return (x[0,4]**2+x[0,5]**2+x[0,6]**2+x[0,7]**2) + 2*(x[0,1]**2+x[0,2]**2+x[0,3]**2) + m1*(x[0,0]**2+x[0,1]**2+x[0,2]**2+x[0,3]**2) + m2*x[0,0]**4 + (x[0,1]**2+x[0,2]**2+x[0,3]**2)*(6*x[0,0]**2+x[0,1]**2+x[0,2]**2+x[0,3])*m2

def toy_model(x, m1=1, m2=50, m3=1):
    '''This is the toy model we have built where we can adjust the size of the interraction term.'''
    return np.matrix([x[0,2], x[0,3], -(m1*2*x[0,0]+m2*2*x[0,0]*x[0,1]**2+m3*4*x[0,0]**3), -(m1*2*x[0,1]+m2*2*x[0,1]*x[0,0]**2+m3*4*x[0,1]**3)])

def hamiltonian_toy_model(x, m1=1, m2=50, m3=1):
    return (x[0,2]**2 + x[0,3]**2) + m1*(x[0,0]**2 + x[0,1]**2) + m2*(x[0,0]**2 * x[0,1]**2) + m3*(x[0,0]**4 + x[0,1]**4)

def initial_condition_toy_model(energy, m1=1, m2=50, m3=1):
    px = 0
    py = 0
    x = np.random.normal()
    polynomial = [m3, 0, m2*x**2 + m1, 0, -energy+px**2+py**2+m3*x**4+m1*x**2]
    solution_list = real_roots(polynomial)
    if len(solution_list) != 0:
        r = solution_list[random.randint(0, len(solution_list)-1)]
        return np.matrix([x, r, px, py, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    else:    
        return initial_condition_toy_model(energy)
    
def new_model(x, m1=1, m2=50, dim=2):
    X = make_matrix(x[0, :dim**2], (dim, dim))
    Y = make_matrix(x[0, dim**2:2*dim**2], (dim, dim))
    PX = make_matrix(x[0, 2*dim**2:3*dim**2], (dim, dim))
    PY = make_matrix(x[0, 3*dim**2:4*dim**2], (dim, dim))

    L = SU2_generators((dim-1)/2)
    L_x, L_y, L_z = L[0], L[1], L[2]

    X_dot = PX.T.flatten()
    Y_dot = PY.T.flatten()
    PX_dot = (2*L_x@X@L_x + 2*L_y@X@L_y + 2*L_z@X@L_z -(m1+dim**2/2-1/2)*X -(m2/2)*(Y@X@X + X@X@Y)).T.flatten()
    PY_dot = (2*L_x@Y@L_x + 2*L_y@Y@L_y + 2*L_z@Y@L_z -(m1+dim**2/2-1/2)*Y -(m2/2)*(X@Y@Y + Y@Y@X)).T.flatten()
    result = np.concatenate((X_dot, Y_dot, PX_dot, PY_dot), axis=1)
    return result

def hamiltonian_new_model(x, m1=1, m2=50, dim=2):
    X = make_matrix(x[0, :dim**2], (dim, dim))
    Y = make_matrix(x[0, dim**2:2*dim**2], (dim, dim))
    PX = make_matrix(x[0, 2*dim**2:3*dim**2], (dim, dim))
    PY = make_matrix(x[0, 3*dim**2:4*dim**2], (dim, dim))

    L = SU2_generators((dim-1)/2)
    L_x, L_y, L_z = L[0], L[1], L[2]

    return 1/2*np.trace(PX@PX + PY@PY + m1*(X@X + Y@Y) + m2*(X@X@Y@Y) 
                        - (commutation(L_x,X))@(commutation(L_x,X)) - (commutation(L_y,X))@(commutation(L_y,X)) - (commutation(L_z,X))@(commutation(L_z,X)) 
                        - (commutation(L_x,Y))@(commutation(L_x,Y)) - (commutation(L_y,Y))@(commutation(L_y,Y)) - (commutation(L_z,Y))@(commutation(L_z,Y)))

def initial_condition_new_model(energy, m1=0.1, m2=0.1, dim=2):
    X = give_random_matrix(dim)
    Y = give_random_matrix(dim)
    Y[0, 0] = np.random.normal()
    PX = zero_matrix((dim,dim))
    PY = zero_matrix((dim,dim))

    L = SU2_generators((dim-1)/2)
    L_x, L_y, L_z = L[0], L[1], L[2]   

    zeroth_order = hamiltonian_new_model(np.concatenate((X.flatten(), Y.flatten(), PX.flatten(), PY.flatten()), axis=1))
    first_order = (m2/2)*(X[0,:]@Y@Y[:,0] + Y[0,:]@Y@X[:,0])[0,0] - 2*(L_x[0,:]@X@L_x[:,0] + L_y[0,:]@X@L_y[:,0] + L_z[0,:]@X@L_z[:,0])[0,0]
    second_order = -L_x[0,0]*L_x[0,0]-L_y[0,0]*L_y[0,0]-L_z[0,0]*L_z[0,0] + (m1/2 + dim**2/4 - 1/4) + (m2/2)*Y[0,0]*Y[0,0]
    polynomial = [second_order.real, first_order.real, zeroth_order.real-energy]
    solution_list = real_roots(polynomial)
    if len(solution_list) != 0:
        r = solution_list[random.randint(0, len(solution_list)-1)]
        X[0, 0] = r
        return split_real_complex(np.concatenate((X.flatten(), Y.flatten(), PX.flatten(), PY.flatten()), axis=1).tolist()[0], dim)
    else:    
        return initial_condition_new_model(energy)
    
def initial_condition_new_model_alternative(energy, m1=1, m2=10, dim=2):
    L = SU2_generators((dim-1)/2)
    L_x, L_y, L_z = L[0], L[1], L[2]
    #alpha = random.uniform(-1/(dim), 1/(dim))
    alpha = np.random.normal(loc=math.sqrt(energy)/dim)
    polynomial = [(m1+2)*(dim*dim*dim-dim)/24+m2*(alpha*alpha)*((dim-1)*(6*(dim-1)*(dim-1)*(dim-1)*(dim-1)+30*(dim-1)*(dim-1)*(dim-1)+40*(dim-1)*(dim-1)-16))/960, 0, alpha*alpha*(m1+2)*(dim*dim*dim-dim)/24-energy]
    solution_list = real_roots(polynomial)
    if len(solution_list) != 0:
        beta = solution_list[random.randint(0, len(solution_list)-1)]
        return [alpha, beta]
    else:
        return initial_condition_new_model_alternative(energy)

def lorenz_model(x, sigma=10, beta=8/3, rho=28):
    x_dot = sigma*(x[0,1] - x[0,0])
    y_dot = x[0,0]*(rho - x[0,2]) - x[0,1]
    z_dot = x[0,0]*x[0,1] - beta*x[0,2]
    return np.matrix([x_dot, y_dot, z_dot])

def random_energy_initial(E, n=2, k=1, m=1):
    pi = math.pi
    r = [np.random.normal(), np.random.normal(), np.random.normal()]

    X = np.sqrt(E)*r[0]/np.sqrt(r[0]**2+r[1]**2+r[2]**2)
    Y = np.sqrt(E)*r[1]/np.sqrt(r[0]**2+r[1]**2+r[2]**2)
    Z = np.sqrt(E)*r[2]/np.sqrt(r[0]**2+r[1]**2+r[2]**2)

    Pq = np.sqrt(4*(n-1)*X**2)

    Pr = np.sqrt(4*(n-1)*Y**2)

    p = [8*pi**2*(n-1)/k**2,0,-4*pi*(n-1)*m/k,0,(n-1)*m**2/2,0,-Z**2]
    p_list = real_roots(p)

    if len(p_list) != 0:
        r = p_list[random.randint(0, len(p_list)-1)]
        result = np.concatenate((np.matrix([0, Pq, r, Pr]), np.matrix(np.identity(4).flatten())), axis=1).tolist()[0]
        return result
    else:
        random_energy_initial(E)


if __name__ == '__main__':
    print(str(initial_condition_new_model_alternative(20)).replace('[','{').replace(']','}'))
