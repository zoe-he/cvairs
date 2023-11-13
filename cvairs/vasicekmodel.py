import csv
import random
import statistics
from math import exp, sqrt


class VasicekModel:
    """Vašíček model (discrete)"""

    def __init__(self, mu, kappa, sigma, r0, N, dt, T, dtau, Tau, seed=777):
        self.mu = mu
        self.kappa = kappa
        self.sigma = sigma
        self.r0 = r0
        self.N = int(N)
        self.dt = dt
        self.T = int(T)
        self.dtau = dtau
        self.Tau = int(Tau)

        self.I = int(T / dt)
        self.t = [i * dt for i in range(self.I + 1)]
        self.J = int(Tau / dtau)
        self.tau = [j * dtau for j in range(self.J + 1)]

        sample_size = self.N * self.I
        random.seed(seed)
        # Standard normal distribution (mu=0.0, sigma=1.0).
        snd = statistics.NormalDist()
        # Calculate samples with inverse CDF.
        self.omega = [snd.inv_cdf(random.random()) for _ in range(sample_size)]

    def to_csv(self, file, rates=None, prices=None):
        if rates is None or prices is None:
            rates = self.short_term_interest_rates()
            prices = self.zero_coupon_bond_prices(rates)
        with open(file, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(['Scen', 'Time', 'Short Rate'] + [f'Price Tau {j}' for j in range(1, self.J + 1)])
            for n in range(self.N):
                for i in range(self.I + 1):
                    writer.writerow([n, self.t[i], rates[n][i]] + [prices[n][i][j] for j in range(self.J)])

    def short_term_interest_rates(self):
        """Generate short-term interest rates.

        r(n, t_i) = κ · μ · ∆t + (1 − κ · ∆t) · r(n, t_(i−1)) + σ · √∆t ω(n, t_i)
        """
        r = []
        for n in range(self.N):
            r.append([self.r0])
            for i in range(1, self.I + 1):
                r[n].append(self.kappa * self.mu * self.dt + (1 - self.kappa * self.dt) * r[n][i - 1] +
                            self.sigma * sqrt(self.dt) *
                            self.omega[self.I * n + (i - 1)])  # Flattened 2D list indexing.
        return r

    def zero_coupon_bond_prices(self, rates=None):
        """Calculate zero-coupon bond prices.

        z(n, t_i, τ_j) = A(τ_j) · e^(−B(τ_j)·r(n,t_i))
        """
        rates = self.short_term_interest_rates() if rates is None else rates
        return [[[self._z(n, i, j, rates)
                  for j in range(1, self.J + 1)]
                 for i in range(self.I + 1)]
                for n in range(self.N)]

    def _z(self, n, i, j, rates):
        btau = self._b(j)
        atau = self._a(j, btau)
        return atau * exp(-btau * rates[n][i])

    def _b(self, j):
        return (1 - exp(-self.kappa * self.tau[j])) / self.kappa

    def _a(self, j, btau):
        kappa_sq = self.kappa ** 2
        sigma_sq = self.sigma ** 2
        return exp((((btau - self.tau[j]) * (kappa_sq * self.mu - sigma_sq / 2)) / kappa_sq) -
                   ((sigma_sq * btau ** 2) / (4 * self.kappa)))
