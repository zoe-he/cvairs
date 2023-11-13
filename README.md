# cvairs

Credit valuation adjustment (CVA) for interest rate swap (IRS):

* Generate short‑term interest rates to calculate zero‑coupon bond prices using the discrete‑time Vašíček model.
* Simulate swap market‑to‑market values estimate expected and potential future exposures of a counterparty entering an IRS.
* Interpolate a credit curve using cubic spline to calculate probabilities of default.
* Calculate a CVA for a counterparty entering an IRS.

## Roadmap

* [x] Reorganize repository as package ([dfc01e4](https://github.com/zoe-he/cvairs/commit/dfc01e4bc313f5c1e3be170bf8ddb548cb269f5e)).
* [x] Add test data and sample script ([27bc5d7](https://github.com/zoe-he/cvairs/commit/27bc5d761417f3dd2cc6b22c502320bc26ef24d7)).
* [x] Add sample Jupyter notebook ([1270df1](https://github.com/zoe-he/cvairs/commit/1270df1267ba2aadc6d82135c31786c72c9f8e74)).
* [ ] Utilize scientific libraries over standard library, e.g., NumPy, SciPy, Matplotlib, and pandas.
* [ ] Add database support via SQLAlchemy with SQLite sample.
