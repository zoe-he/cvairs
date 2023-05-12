#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['CCR', 'CVA', 'IRSExposure', 'VasicekModel']

import csv
import statistics  # NormalDist and fmean require Python 3.8+
from math import exp, fsum, sqrt
from random import random


class VasicekModel:
    """Vašíček model (discrete)"""

    def __init__(self, mu, kappa, sigma, r0, N, dt, T, dtau, Tau):
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
        snd = statistics.NormalDist()  # Standard normal distribution (mu=0.0, sigma=1.0).
        self.omega = [snd.inv_cdf(random()) for _ in range(sample_size)]  # Calculate samples with inverse CDF.

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


class IRSExposure:
    """Interest rate swap (IRS) exposure (receive fixed by default)"""

    def __init__(self, s, k, vasicek_model):
        self.s = s
        self.k = k

        self.N = vasicek_model.N
        self.dt = vasicek_model.dt
        self.I = vasicek_model.I
        self.t = vasicek_model.t
        self.z = vasicek_model.zero_coupon_bond_prices()

    def to_csv(self, file, mtms=None, *, swaption_type='receiver'):
        mtms = self.market_to_market_simulations(swaption_type=swaption_type) if mtms is None else mtms
        ee = self.expected_exposures(mtms)
        pfe = self.potential_future_exposures(mtms)
        with open(file, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(['Time', 'EE', f'PFE({self.k})'] + [f'Scen {n}' for n in range(1, self.N + 1)])
            for i in range(self.I + 1):
                writer.writerow([self.t[i], ee[i], pfe[i]] + mtms[i])

    def market_to_market_simulations(self, *, swaption_type='receiver'):
        """Simulate swap market-to-market (MTM) values."""
        return [[self._swap_mtm(n, i, swaption_type) for n in range(self.N)] for i in range(self.I + 1)]

    def expected_exposures(self, mtms):
        """Calculate expected exposures (EE)."""
        return [self._ee(i, mtms) for i in range(self.I + 1)]

    def _ee(self, i, mtms):
        return fsum([mtm for mtm in mtms[i] if mtm > 0]) / self.N

    def potential_future_exposures(self, mtms):
        """Calculate potential future exposures (PFE)."""
        return [self._pfe(i, mtms) for i in range(self.I + 1)]

    def _pfe(self, i, mtms):
        kth_index = int(self.k * self.N) - 1  # Trunc floats, like Excel's SMALL function.
        sorted_mtms = sorted(mtms[i])
        return sorted_mtms[kth_index]

    def _fixed_leg(self, n, i):
        return self.s * self.dt * fsum(self.z[n][i][:self.I - i])

    def _floating_leg(self, n, i):
        return 1 - self.z[n][i][self.I - i - 1] if i != self.I else 0.0

    def _swap_mtm(self, n, i, swaption_type):
        if swaption_type == 'receiver':
            return self._fixed_leg(n, i) - self._floating_leg(n, i)
        if swaption_type == 'payer':
            return self._floating_leg(n, i) - self._fixed_leg(n, i)
        raise ValueError(f'swaption_type={swaption_type} is not allowed.')


class CCR:
    """Counterparty credit risk (CCR)"""

    @classmethod
    def to_csv(cls, file, credit_curve, dt, LGD):
        x, y = zip(*credit_curve)
        t, y_interp = cls.cubic_spline(x, y, dt)
        pds = cls.default_probabilities(t, y_interp, LGD)
        with open(file, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(['t', 'Credit Spread (t)', 'Default Prob'])
            for i in range(len(t)):
                writer.writerow([t[i], y_interp[i], pds[i]])

    @classmethod
    def cubic_spline(cls, x, y, dt):
        """Calculate cubic spline."""
        T = x[-1]
        I = int(T / dt)
        t = [i * dt for i in range(I + 1)]
        y2 = cls._spline(x, y)
        y_interp = [cls._interp(x, y, y2, xv) for xv in t]
        return t, y_interp

    @staticmethod
    def _spline(xl, yl, yp1=1.e99, ypn=1.e99):
        # Translated from Numerical Recipes 3E (C++).
        # interp_1d.h: void Spline_interp::sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn)
        n = len(xl)
        assert n == len(yl)
        y2 = [0.0] * n
        u = [0.0] * (n - 1)
        if yp1 > 0.99e99:  # Lower boundary condition type is natural.
            y2[0] = 0.0
            u[0] = 0.0
        else:  # Otherwise, has specified first derivative value.
            y2[0] = -0.5
            u[0] = (3.0 / (xl[1] - xl[0])) * ((yl[1] - yl[0]) / (xl[1] - xl[0]) - yp1)
        for i in range(1, n - 1):
            sig = (xl[i] - xl[i - 1]) / (xl[i + 1] - xl[i - 1])
            p = sig * y2[i - 1] + 2.0
            y2[i] = (sig - 1.0) / p
            u[i] = (yl[i + 1] - yl[i]) / (xl[i + 1] - xl[i]) - (yl[i] - yl[i - 1]) / (xl[i] - xl[i - 1])
            u[i] = (6.0 * u[i] / (xl[i + 1] - xl[i - 1]) - sig * u[i - 1]) / p
        if ypn > 0.99e99:  # Upper boundary condition type is natural.
            qn = 0.0
            un = 0.0
        else:  # Otherwise, has specified first derivative value.
            qn = 0.5
            un = (3.0 / (xl[n - 1] - xl[n - 2])) * (ypn - (yl[n - 1] - yl[n - 2]) / (xl[n - 1] - xl[n - 2]))
        y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)
        for k in range(n - 2, -1, -1):
            y2[k] = y2[k] * y2[k + 1] + u[k]
        return y2

    @staticmethod
    def _interp(xl, yl, y2l, x):
        # Translated from Numerical Recipes 3E (C++).
        # interp_1d.h: Doub Spline_interp::rawinterp(Int jl, Doub x), Int Base_interp::locate(const Doub x)
        n = len(xl)
        if n < 2:
            raise ValueError('xl length must be greater than or equal to 2.')
        klo = 0
        khi = n - 1
        while khi - klo > 1:
            k = (khi + klo) >> 1
            if xl[k] > x:
                khi = k
            else:
                klo = k
        h = xl[khi] - xl[klo]
        if h == 0.0:
            raise ValueError('xl values must be distinct.')
        a = (xl[khi] - x) / h
        b = (x - xl[klo]) / h
        return a * yl[klo] + b * yl[khi] + ((a ** 3 - a) * y2l[klo] + (b ** 3 - b) * y2l[khi]) * (h ** 2) / 6.0

    @classmethod
    def default_probabilities(cls, times, interpolated_credit_spreads, LGD):
        """Calculate probabilities of default (PD)."""
        pds = [0.0]
        pds += [cls._pd(times, interpolated_credit_spreads, LGD, i) for i in range(1, len(times))]
        return pds

    @staticmethod
    def _pd(t, S, LGD, i):
        return exp(-S[i - 1] * t[i - 1] / 10_000 / LGD) - exp(-S[i] * t[i] / 10_000 / LGD)


class CVA:
    """Credit value adjustment (CVA)"""

    @staticmethod
    def cva(LGD, EE, PD):
        return -LGD * fsum(ee * pd for ee, pd in zip(EE, PD))

    @staticmethod
    def risky_duration(S, LGD):
        h = S[-1] / 10_000 / LGD
        return (1 - exp(-h * 10)) / h

    @staticmethod
    def cva_bps_approx1(CVA, RD):
        return CVA * 10_000 / RD

    @staticmethod
    def epe(EE):
        return statistics.fmean(EE)

    @staticmethod
    def cva_bps_approx2(EPE, S):
        return -EPE * S[-1]

    @staticmethod
    def to_csv(file, time, interpolated_credit_spread, probabilities_of_default, expected_exposures, cva_values):
        with open(file, 'w', newline='') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(['Time (years)', 'Interpolated Credit Curve (Spread in bps)', 'Default Prob', 'EE',
                             'Name', 'Value'])
            n = len(cva_values)
            for i in range(len(time)):
                row = [time[i], interpolated_credit_spread[i], probabilities_of_default[i], expected_exposures[i]]
                if i < n:
                    row.extend(cva_values[i])
                writer.writerow(row)


def _params_from_csv(file):
    with open(file, newline='') as params_file:
        reader = csv.reader(params_file)
        next(reader)  # Skip main params header.
        params = []  # [Value: float]
        for row in reader:
            if not row[1]:
                break
            params.append(float(row[1]))
        next(reader)  # Skip credit curve header.
        credit_curve = [(float(row[0]), float(row[1])) for row in reader]  # [(Time: years, Credit spread: bps)]
        return params, credit_curve
