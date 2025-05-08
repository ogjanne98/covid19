from datetime import date, datetime
import numpy as np
import matplotlib.pyplot as plt

class Beta:
    """Class for creating piecewise constant functon from pre-defined
    date-value data"""
    def __init__(self, filename):
        self.B = []
        self.dates = [0]
        with open(filename, "r") as infile:
            for _ in range(2):
                infile.readline()
            for line in infile:
                data = line.split()
                self.B.append(float(data[2]))
                day, month, year = [int(n) for n in data[0].split(".")]
                date0 = date(year, month, day)
                day, month, year = [int(n) for n in data[1].split(".")]
                date1 = date(year, month, day)
                delta = date1-date0
                n_days = delta.days
                self.dates.append(self.dates[-1] + int(n_days))

    def __call__(self, T):
        """T represents number of days starting from 15.02.2020"""
        for k in range(len(self.dates) - 1):
            if (self.dates[k] <= T) and (T <= self.dates[k+1]):
                return self.B[k]
        return self.B[-1]

    def plot(self, T):
        self.t = np.linspace(0,T, 1001)
        self.y = np.zeros(len(self.t))
        index = 0
        for k in self.t:
            self.y[index] = self(k)
            index = index + 1
        plt.plot(self.t,self.y)
        plt.savefig("beta.pdf")
        plt.show()

if __name__ == "__main__":
    from outbreak import outbreak_Norway
    beta = Beta("beta_values.txt")
    beta.plot(1000)
    outbreak_Norway(beta,1000,1)
