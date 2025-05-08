from SEIR import *
import numpy as np
import matplotlib.pyplot as plt
from Beta import Beta

class SEIRimport(ProblemSEIR):
    """Add influx of exposed people, e.g. from tourists and travelers, by parameter sgma"""
    def __init__(self, beta, r_ia = 0.1, r_e2=1.25, lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2, sgma=10):
        ProblemSEIR.__init__(self, beta, r_ia, r_e2, lmbda_1, lmbda_2, p_a, mu)
        self.sgma = sgma

    def __call__(self, t, y):
        """Add influx to E2 category"""
        dP = ProblemSEIR.__call__(self, t, y)
        dP[2] = dP[2] + self.sgma
        return  dP

def outbreak_Norway(beta, n_days, dt):
    S_0 = 5e6; E2_0 = 100
    problem = SEIRimport(beta)
    problem.set_initial_condition(S_0, E2_0)
    solver = SolverSEIR(problem, n_days, dt)
    solver.solve()
    solver.plot(["I", "Ia"])

if __name__ == "__main__":
    beta = Beta("beta_values.txt")
    beta.plot(1000)
    outbreak_Norway(beta,1000,1)
