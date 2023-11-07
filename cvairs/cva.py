import csv
import statistics
from math import exp, fsum


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
