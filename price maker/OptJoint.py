import numpy as np
import gurobipy as gp
from gurobipy import GRB

def jointOpt_taker(T, p, w_true, s_ub, s_lb, s_init, eta_c, eta_d):
    m = gp.Model("Joint")
    w = m.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output")
    d = m.addMVar(shape= T, name="storage discharge")
    c = m.addMVar(shape= T, name="storage charge")
    s = m.addMVar(shape= T, name="storage SOC")
    o = m.addMVar(shape= T, lb=-GRB.INFINITY, name="joint operation")
    x = m.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    m.addConstr(o == x + w, name="joint supply")
    m.addConstr(d-c >= x + w - w_true, name="joint storagesupply")

    m. setObjective(p @ o, GRB.MAXIMIZE)
    m.addConstr(s <= s_ub, name="soc_ub")
    m.addConstr(s >= s_lb, name="soc_lb")
    add_mat = np.eye(T, k=-1)
    # print(add_mat)
    m.addConstr(s == add_mat @ s - (1/eta_d)*d + eta_c* c + s_init)
    m.optimize()
    return w.X, x.X, s.X, m.ObjVal

def jointOpt_taker_noArb(T, p, w_true, s_ub, s_lb, s_init, eta_c, eta_d):
    m_no = gp.Model("Joint_no_arbitrage")
    w3 = m_no.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output")
    d3 = m_no.addMVar(shape= T, name="storage discharge")
    c3 = m_no.addMVar(shape= T, name="storage charge")
    x3 = m_no.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    s3 = m_no.addMVar(shape= T, name="storage SOC")
    o3 = m_no.addMVar(shape= T, lb=-GRB.INFINITY, name="joint operation")
    m_no.addConstr(d3 - c3 >= x3 + w3 - w_true, name="joint storage")
    m_no.addConstr(o3 == x3 + w3, name="joint supply")
    ## no arbitrage
    m_no.addConstr(o3 >= np.zeros(T), name="no buying")
    # m_no.addConstr(c3 <= w3, name="no buying")
    ########
    m_no. setObjective(p @ o3, GRB.MAXIMIZE)
    m_no.addConstr(s3 <= s_ub, name="soc_ub")
    m_no.addConstr(s3 >= s_lb, name="soc_lb")
    add_mat = np.eye(T, k=-1)
    m_no.addConstr(s3 == add_mat @ s3 - (1/eta_d)*d3 + eta_c* c3 + s_init)
    m_no.optimize()
    return w3.X, x3.X, s3.X, m_no.ObjVal


def jointOpt_linear(T, alpha, beta, w_true, s_ub, s_lb, s_init, eta_c, eta_d):
    m = gp.Model("Joint")
    m.Params.LogToConsole = 0
    w = m.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output")
    d = m.addMVar(shape= T, name="storage discharge")
    c = m.addMVar(shape= T, name="storage charge")
    s = m.addMVar(shape= T, name="storage SOC")
    o = m.addMVar(shape= T, lb=-GRB.INFINITY, name="joint operation")
    x = m.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    m.addConstr(o == x + w, name="joint supply")
    m.addConstr(d-c >= x + w - w_true, name="joint storage")

    m. setObjective(alpha @ o - o @ beta @ o, GRB.MAXIMIZE)
    m.addConstr(s <= s_ub, name="soc_ub")
    m.addConstr(s >= s_lb, name="soc_lb")
    add_mat = np.eye(T, k=-1)
    # print(add_mat)
    m.addConstr(s == add_mat @ s - (1/eta_d)*d + eta_c* c + s_init)
    m.optimize()
    return w.X, x.X, s.X, m.ObjVal

def jointOpt_linear_noArb(T, alpha, beta, w_true, s_ub, s_lb, s_init, eta_c, eta_d):
    m = gp.Model("Joint_noArb")
    m.Params.LogToConsole = 0
    w = m.addMVar(shape= T, lb=-GRB.INFINITY, name="wind output")
    d = m.addMVar(shape= T, name="storage discharge")
    c = m.addMVar(shape= T, name="storage charge")
    s = m.addMVar(shape= T, name="storage SOC")
    o = m.addMVar(shape= T, lb=-GRB.INFINITY, name="joint operation")
    x = m.addMVar(shape= T, lb=-GRB.INFINITY, name="storage operation")
    m.addConstr(o == x + w, name="joint supply")
    m.addConstr(d-c >= x + w - w_true, name="joint storage")
    ## no arbitrage
    m.addConstr(o >= np.zeros(T), name="no buying")

    m. setObjective(alpha @ o - o @ beta @ o, GRB.MAXIMIZE)
    m.addConstr(s <= s_ub, name="soc_ub")
    m.addConstr(s >= s_lb, name="soc_lb")
    add_mat = np.eye(T, k=-1)
    # print(add_mat)
    m.addConstr(s == add_mat @ s - (1/eta_d)*d + eta_c* c + s_init)
    m.optimize()
    return w.X, x.X, s.X, m.ObjVal