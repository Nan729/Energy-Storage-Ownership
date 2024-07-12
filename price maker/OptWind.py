import numpy as np
import gurobipy as gp
from gurobipy import GRB

def windOpt_taker(T, p, w_true):
    m2_w =  gp.Model("DisJoint-wind")
    w2 = m2_w.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output 2")
    m2_w.setObjective(p @ w2, GRB.MAXIMIZE)
    m2_w.addConstr( w2 <= w_true)
    m2_w.optimize()
    return w2.X, m2_w.ObjVal

def windOpt_linear(T, alpha, beta, w_true):
    m2_w =  gp.Model("DisJoint-wind")
    m2_w.Params.LogToConsole = 0
    w2 = m2_w.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output 2")
    m2_w.setObjective(alpha @ w2 - w2 @ beta @ w2, GRB.MAXIMIZE)
    m2_w.addConstr( w2 <= w_true)
    m2_w.optimize()
    return w2.X, m2_w.ObjVal