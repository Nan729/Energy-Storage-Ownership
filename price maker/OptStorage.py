import numpy as np
import gurobipy as gp
from gurobipy import GRB

def storageOpt_taker(T, p, s_ub, s_lb, eta_d, eta_c, s_init):
    m2_s =  gp.Model("DisJoint-storage")
    d2 = m2_s.addMVar(shape= T, name="storage discharge 2")
    c2 = m2_s.addMVar(shape= T, name="storage charge 2")
    s2 = m2_s.addMVar(shape= T, name="storage SOC 2")
    x2 = m2_s.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    m2_s.addConstr( x2 == d2 - c2)
    m2_s.setObjective(p @ x2, GRB.MAXIMIZE)
    m2_s.addConstr( s2 <= s_ub)
    m2_s.addConstr( s2 >= s_lb)
    add_mat = np.eye(T, k=-1)
    m2_s.addConstr( s2 == add_mat @ s2 - (1/eta_d)*d2 + eta_c* c2 + s_init)
    m2_s.optimize()
    return x2.X, m2_s.ObjVal

def storageOpt_linear(T, alpha, beta, s_ub, s_lb, eta_d, eta_c, s_init):
    m2_s =  gp.Model("DisJoint-storage")
    m2_s.Params.LogToConsole = 0
    d2 = m2_s.addMVar(shape= T, name="storage discharge 2")
    c2 = m2_s.addMVar(shape= T, name="storage charge 2")
    s2 = m2_s.addMVar(shape= T, name="storage SOC 2")
    x2 = m2_s.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    m2_s.addConstr( x2 == d2 - c2)
    m2_s.setObjective(alpha @ x2- x2@ beta @x2, GRB.MAXIMIZE)
    m2_s.addConstr( s2 <= s_ub)
    m2_s.addConstr( s2 >= s_lb)
    add_mat = np.eye(T, k=-1)
    m2_s.addConstr( s2 == add_mat @ s2 - (1/eta_d)*d2 + eta_c* c2 + s_init)
    m2_s.optimize()
    return x2.X, m2_s.ObjVal