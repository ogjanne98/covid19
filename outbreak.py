from SEIR import *

def outbreak_Norway(beta, n_days, dt):
    S_0 = 5e6; E2_0 = 100
    problem = ProblemSEIR(beta)
    problem.set_initial_condition(S_0, E2_0)
    solver = SolverSEIR(problem, n_days, dt)
    solver.solve()
    solver.plot(["I", "Ia"])
    max_infected = max(solver.y[:,3])
    print(max_infected*0.05)

def B(t):
    if t < 30:
        return 0.33
    else:
        return 0.083


outbreak_Norway(0.33, 150, 1)
outbreak_Norway(B, 150, 1)
