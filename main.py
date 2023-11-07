import csv

from cvairs import CCR, CVA, IRSExposure, VasicekModel


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


if __name__ == '__main__':
    params, credit_curve = _params_from_csv('data/params.csv')
    dt, LGD = params[5], params[11]

    model = VasicekModel(*params[:9])

    exposure = IRSExposure(*params[9:11], model)
    mtms = exposure.market_to_market_simulations()
    EE = exposure.expected_exposures(mtms)

    x, y = zip(*credit_curve)
    t, S = CCR.cubic_spline(x, y, dt)
    PD = CCR.default_probabilities(t, S, LGD)

    cva = CVA.cva(LGD, EE, PD)
    risky_duration = CVA.risky_duration(S, LGD)
    cva_bps_approx1 = CVA.cva_bps_approx1(cva, risky_duration)
    epe = CVA.epe(EE)
    cva_bps_approx2 = CVA.cva_bps_approx2(epe, S)
    cva_values = [('CVA', cva), ('Risky Duration', risky_duration), ('CVA (bps) Approximation 1', cva_bps_approx1),
                  ('EPE', epe), ('CVA (bps) Approximation 2', cva_bps_approx2)]
    CVA.to_csv('out.csv', t, S, PD, EE, cva_values)
