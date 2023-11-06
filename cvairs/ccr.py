import csv
from math import exp


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
