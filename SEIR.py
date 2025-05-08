import numpy as np
import matplotlib.pyplot as plt
from ODEsolver import *

class ProblemSEIR:
    """Implement the SEEIIR model of infectious disease with callable beta(t)"""
    def __init__(self, beta, r_ia = 0.1, r_e2=1.25, \
    lmbda_1=0.33, lmbda_2=0.5, p_a=0.4, mu=0.2):
        if isinstance(beta, (float,int)): # check if beta constant or callable
            self.beta = lambda t: beta
        elif callable(beta):
            self.beta = beta
        else:
            msg = "Parameter \'beta\' must be a constant or a function of time"
            raise ImplementationError(msg)
        self.r_ia = r_ia
        self.r_e2 = r_e2
        self.lmbda_1 = lmbda_1
        self.lmbda_2 = lmbda_2
        self.p_a = p_a
        self.mu = mu

    def set_initial_condition(self, S_0, E2_0):
        self.initial_condition = [S_0, 0, E2_0, 0, 0, 0]

    def get_population(self):
        try:
            return sum(self.initial_condition)
        except:
            msg = "Set initial condition before calling this method"
            raise NoInitialConditionError(msg)

    def __call__(self, t, y):
        """Make instances callable and define ODE equations"""
        beta = self.beta
        r_ia = self.r_ia
        r_e2 = self.r_e2
        lmbda_1 = self.lmbda_1
        lmbda_2 = self.lmbda_2
        p_a = self.p_a
        mu = self.mu
        S, E1, E2, I, Ia, R = y
        N = self.get_population()
        dS = -beta(t) * S * I / N - r_ia * beta(t) * S * Ia / N \
        - r_e2 * beta(t) * S * E2 / N
        dE1 = beta(t) * S * I / N + r_ia * beta(t) * S * Ia / N \
        + r_e2 * beta(t) * S * E2 / N - lmbda_1 * E1
        dE2 = lmbda_1 * (1 - p_a) * E1 - lmbda_2 * E2
        dI = lmbda_2 * E2 - mu * I
        dIa = lmbda_1 * p_a * E1 - mu * Ia
        dR = mu * (I + Ia)
        return [dS, dE1, dE2, dI, dIa, dR]

class SolverSEIR:
    """Class for solving and plotting pre-defined SEEIIR-problem"""
    def __init__(self, problem, T, dt):
        self.problem = lambda t,y: np.asarray(problem(t,y))
        self.T = T # number of days
        self.dt = dt
        self.y0 = problem.initial_condition

    def solve(self, method=ForwardEuler):
        """Solve with method of choice from ODEsolver class"""
        solver = method(self.problem)
        solver.set_initial_condition(self.y0)
        N = int(self.T//self.dt)
        self.t, self.y = solver.solve([0,self.T], N)

    def plot(self, states):
        """Plot the given categories in states(list of category strings)"""
        states_0 = ["S", "E1", "E2", "I", "Ia", "R"]
        for state in states:
            index = states_0.index(state)
            plt.plot(self.t, self.y[:,index], label=state)
        plt.legend()
        plt.show()

def test_ProblemSEIR():
    instance = ProblemSEIR(beta=0.4, r_ia=0.1, r_e2=1.25, lmbda_1=0.33,
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
    S_0 = 100; E2_0 = 1000
    instance.set_initial_condition(S_0, E2_0)
    assert [S_0, 0, E2_0, 0, 0, 0] == instance.initial_condition
    assert instance.get_population() == 1100

if __name__ == "__main__":
    test_ProblemSEIR()
    S_0 = 5e6
    E2_0 = 100
    problem = ProblemSEIR(beta=0.4)
    problem.set_initial_condition(S_0,E2_0)
    solver = SolverSEIR(problem,T=150,dt=1.0)
    solver.solve()
    solver.plot(["S", "I", "Ia", "R"])
