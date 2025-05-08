import numpy as np
import matplotlib.pyplot as plt
from SEEIIR import *
from ODEsolver import *

# Simulate spread of virus in society, plot results.
def COV19(beta=0.5, r_ia=0.1, r_e2=1.25, lmbda_1=0.33, lmbda_2=1, p_a=0, mu=0.2, d=0.5):
    seeiir = SEEIIR(beta, r_ia, r_e2, lmbda_1, lmbda_2, p_a, mu, d)
    solver = Midpoint(seeiir)
    solver.set_initial_condition([5*(10**6), 0, 100, 0, 0, 0,0])
    t,y = solver.solve([0,365], 365)
    plt.plot(t, y[:,0], label = "S")
    plt.plot(t,y[:,3], label="I")
    plt.plot(t,y[:,5], label="R")
    plt.plot(t,y[:,6], label="D")
    plt.title("Covid-19 spread in Norway")
    plt.legend()
    plt.show()

COV19()    
