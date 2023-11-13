# cvairs

Credit valuation adjustment (CVA) for interest rate swap (IRS):

* Generate short‑term interest rates to calculate zero‑coupon bond prices using the discrete‑time Vašíček model.
* Simulate swap market‑to‑market values estimate expected and potential future exposures of a counterparty entering an IRS.
* Interpolate a credit curve using cubic spline to calculate probabilities of default.
* Calculate a CVA for a counterparty entering an IRS.

## Roadmap

* [x] Reorganize repository as package.
* [x] Add test data and sample script.
* [x] Add sample Jupyter notebook.
* [ ] Utilize scientific libraries over standard library, e.g., NumPy, SciPy, Matplotlib, and pandas.
* [ ] Add database support via SQLAlchemy with SQLite sample.
