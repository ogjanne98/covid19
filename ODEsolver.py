import numpy as np
import matplotlib.pyplot as plt

class ODEsolver: 
    """Class for solving ordinary differential equations numerically, superclass
    not to be instansiated"""
    def __init__(self, f):
        self.f = lambda t,y: np.asarray(f(t,y)) 

    def set_initial_condition(self, y0):
        """Set initial condition and determine number of equations"""
        if isinstance(y0, float) or isinstance(y0, int):
            self.y0 = float(y0)
            self.neq = 1 # number of equations
        else:
            self.y0 = np.asarray(y0)
            self.neq = len(y0)

    def advance(self):
        raise NotImplementedError(\
        "No advance method implemented in superclass")


    def solve(self, interval, N):
        """Solve ODE numerically, interval = [t0,t1], N number of time steps"""
        t0, t1 = interval
        self.dt = (t1-t0)/N
        self.t = np.linspace(t0, t1, N+1)
        if self.neq == 1:
            self.y = np.zeros(N+1)
        else:
            self.y = np.zeros((N+1, self.neq))
        self.y[0] = self.y0
        for n in range(N):
            self.n = n
            self.y[n+1] = self.advance()
        return self.t, self.y

class ForwardEuler(ODEsolver):
    def advance(self):
        f, y, dt, t = self.f, self.y[self.n], self.dt, self.t[self.n]
        return y + dt*f(t, y)

class Midpoint(ODEsolver):
    def advance(self):
        f, y, dt, t = self.f, self.y[self.n], self.dt, self.t[self.n]
        k1 = f(t,y)
        k2 = f(t + dt/2, y + k1*(dt/2))
        return y + dt*k2

class Trapezoidal(ODEsolver):
    def advance(self):
        f, y, dt, t = self.f, self.y[self.n], self.dt, self.t[self.n]
        k1 = f(t,y)
        k2 = f(t + dt, y + dt*k1)
        return y + (dt/2)*(k1 + k2)
