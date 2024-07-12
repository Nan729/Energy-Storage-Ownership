import numpy as np
import gurobipy as gp
from gurobipy import GRB

from OptWind import windOpt_linear
from OptStorage import storageOpt_linear

def disjointOpt_linear(T, alpha, beta, w_true, s_ub, s_lb, s_init, eta_c, eta_d, eps, kUb):
    #### How to compute the convergence
    e_diff = np.ones(T)*10
    # eps = 1e-2
    ## 0. initialize 
    w_sol, _ = windOpt_linear(T, alpha, beta, w_true)
    x_sol = np.zeros(T)
    ## loop: wind-> storage -> wind ->...-> converge
    k = 0
    while np.linalg.norm(e_diff)>eps and k< kUb:
        w_prev = w_sol
        x_prev = x_sol
        alpha_new = alpha-w_sol@beta
        x_sol, m_disjoint_storage = storageOpt_linear(T, alpha_new, beta, s_ub, s_lb, eta_d, eta_c, s_init)
        alpha_new = alpha-x_sol@beta
        w_sol, m_disjoint_wind = windOpt_linear(T, alpha_new, beta, w_true)
        e_diff = np.absolute(x_sol-x_prev) + np.absolute(w_sol-w_prev)
        k += 1 
        # print(k, x_sol, w_sol, e_diff )
    return w_sol, x_sol, m_disjoint_wind, m_disjoint_storage