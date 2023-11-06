import csv
from math import fsum


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
