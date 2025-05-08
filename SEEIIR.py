
class SEEIIR:
    def __init__(self, beta=0.33, r_ia=0.1,
                r_e2=1.25, lmbda_1=0.33,
                lmbda_2=0.5, p_a=0.4, mu=0.2, d= 0.01):
            self.beta = beta
            self.r_ia = r_ia
            self.r_e2 = r_e2
            self.lmbda_1 = lmbda_1
            self.lmbda_2 = lmbda_2
            self.p_a = p_a
            self.mu = mu
            self.d = d 

    def __call__(self,t,y):
            beta = self.beta
            r_ia = self.r_ia
            r_e2 = self.r_e2
            lmbda_1 = self.lmbda_1
            lmbda_2 = self.lmbda_2
            p_a = self.p_a
            mu = self.mu
            d = self.d
            S, E1, E2, I, Ia, R, D = y
            N = sum(y)
            dS = -beta * S * I / N - r_ia * beta * S * Ia / N \
            - r_e2 * beta * S * E2 / N
            dE1 = beta * S * I / N + r_ia * beta * S * Ia / N \
            + r_e2 * beta * S * E2 / N - lmbda_1 * E1
            dE2 = lmbda_1 * (1 - p_a) * E1 - lmbda_2 * E2
            dI = lmbda_2 * E2 - mu * I
            dIa = lmbda_1 * p_a * E1 - mu * Ia
            dR = mu * ((1-d)*I + Ia)
            dD = mu*d*I
            return [dS, dE1, dE2, dI, dIa, dR, dD]
