Below is a cleaner, more structured, and clearer version of your README section. It improves readability, fixes formatting issues, and adds short descriptions for each script.

---

## Setup

This project uses **`renv`** for reproducible package management.

### 1. Install and restore the environment

In R, run:

```r
install.packages("renv")
renv::restore()
```

This will install all required packages with the exact versions used in this project.

---

## How to run the code

The project consists of three independent scripts. Each script can be run separately depending on the task.

### 1. Monte Carlo Simulation – MLE

```r
source("MLE.R")
```

**Description:**
Simulates credit portfolios under a one-factor Gaussian model and estimates asset correlation using **Maximum Likelihood Estimation (MLE)**.
Also analyzes the distribution of the MLE estimators across simulations.

---

### 2. Monte Carlo Simulation – Method of Moments (MoM)

```r
source("MoM.R")
```

**Description:**
Simulates credit portfolios and estimates asset correlation using the **Method of Moments (MoM)**.
Includes several **bias-adjusted estimators** and compares their performance.

---

### 3. Empirical Analysis – IBRD Data

```r
source("Empirical_IBRD.R")
```

**Description:**
Applies both **MLE** and **MoM estimators** to real-world default data (IBRD dataset).
Also includes:

* Comparison with estimators from the `AssetCorr` package
* Basel-style correlation calibration
* Capital (unexpected loss) calculations under different LGD assumptions

---

## Notes

* The scripts are **self-contained** and do not need to be run in a specific order.
* The file `IBRD.xlsx` must be located in the working directory for the empirical analysis.
* Parameter values (e.g., number of simulations, portfolio size) can be adjusted directly within each script.

---

