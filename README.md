# Multi‑Peak Gaussian Fitting with AIC Model Selection

<img width="640" height="480" alt="curve_fit_ex" src="https://github.com/user-attachments/assets/ba4ae486-addf-47bf-ae0b-8ea0c0fa6d3f" />

This project performs multi‑peak Gaussian fitting on noisy data using SciPy’s nonlinear least‑squares optimizer. It detects peaks, builds initial parameter guesses, fits one or more Gaussian components, and compares model complexity using the Akaike Information Criterion (AIC).

My current task focuses on generalizing AIC testing to multi‑peak models, allowing automatic selection of the optimal number of Gaussian peaks.

## Features
- Gaussian and multi‑Gaussian model definitions  
- Synthetic noisy dataset generation  
- Peak detection using `scipy.signal.find_peaks`  
- Automatic parameter initialization from detected peaks  
- Nonlinear curve fitting with `scipy.optimize.curve_fit`  
- AIC calculation for model comparison  
- Plotting of noisy data and fitted curves  

## Requirements
Install dependencies with:

```bash
pip install numpy scipy matplotlib
