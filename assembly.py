
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import minimize
import matplotlib.plot as pyplot

N = 1000
E = 1e9
S = 0.01
F_d = 50
L = 10.0
h = L / (N - 1)
k = E * S / h
g = 100 # penalty

# direct elimination:
def direct_elimination(k, F_d, N, u, f):
    
    a = - np.ones(N - 2)
    b = 2 * np.ones(N - 1)
    b[-1] = 1

    K_eliminated = k * diags([a, b, a], offsets=[-1, 0, 1], format='csr')    

    f_eliminated = np.zeros(N - 1)
    f_eliminated[-1] = F_d

    u_eliminated = spsolve(K_eliminated, f_eliminated)
    f_1 = - k * u_eliminated[0]
    
    u = np.concatenate([[0], u_eliminated])    
    f = np.concatenate([[f_1], f_eliminated])
    
    return u, f

# penalty method:
def penalty_method(k, F_d, N, g, u, f):
    
    a = -np.ones(N - 1)
    b = 2 * np.ones(N)
    b[0] = g + 1
    b[-1] = 1
    
    K = k * diags([a, b, a], offsets=[-1, 0, 1], format='csr')
    
    def objective(x):
        
        f_1 = x[0]
        u = x[1:]
        
        residual = K.dot(u) - np.concatenate([[f_1], np.zeros[N - 2], [F_d]])
        
        penalty = g * f_1**2
        
        return np.linalg.norm(residual)**2 + penalty
    
    x0 = np.zeros(N + 1)
    result = minimize(objective, x0, method='BFGS')
    
    x_opt = result.x
    u = x_opt[1:]
    f = np.concatenate([[x_opt[0]], np.zeros[N - 2], [F_d]])
    
    return u, f

# lagrangian multiplier method:
def lagrangian_multiplier(k, F_d, N, lam, u, f):    
    
    
    
# plot results:
def plot_compare(E, S, L, F_d, N):