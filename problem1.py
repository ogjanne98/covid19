import numpy as np
from SEEIIR import SEEIIR
from ODEsolver import ForwardEuler
import matplotlib.pyplot as plt

def test_SEEIIR():
    instance = SEEIIR(beta=0.4, r_ia=0.1, r_e2=1.25, lmbda_1=0.33,
                       lmbda_2=0.5, p_a=0.4, mu=0.2)
    t = 0
    u = [1,1,1,1,1,1]
    exp = [-0.156666666666, -0.1733333333333, -0.302, 0.3, -0.068, 0.4]
    comp = instance(t,u)
    tol = 1e-10
    msg = "The __call__ function of class SEEIIR returns incorrect value"
    for ex,com in zip(exp,comp):
        bool = abs(ex-com) < tol
        assert bool, msg

def solve_SEEIIR(T,dt,S_0,E2_0):
    y0 = [S_0, 0, E2_0, 0, 0, 0]
    range = [0,T]
    N = T//dt
    f = SEEIIR(beta=0.4, r_ia=0.1, r_e2=1.25,
               lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2)
    solver = ForwardEuler(f)
    solver.set_initial_condition(y0)
    t,y = solver.solve(range,N)
    return t,y

def plot_SEEIIR():
    T=150; dt=1; S_0=5e6; E2_0=100
    t,y = solve_SEEIIR(T,dt,S_0,E2_0)
    groups = {"S": y[:,0], "I": y[:,3], "Ia": y[:,4], "R": y[:,5]}
    for group in groups:
        plt.plot(t, groups[group], label=f"{group}")
    plt.legend()
    plt.show()

plot_SEEIIR()
