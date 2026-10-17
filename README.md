# Multi‑Peak Gaussian Fitting with BIC Model Selection

<img width="640" height="480" alt="curve_fit_ex" src="https://github.com/user-attachments/assets/ba4ae486-addf-47bf-ae0b-8ea0c0fa6d3f" />

This project performs multi‑peak Gaussian fitting on noisy data using SciPy’s nonlinear least‑squares optimizer. It detects peaks, builds initial parameter guesses, fits one or more Gaussian components, and compares multi-peak error using BIC (Bayesian Information Criterion).

## Features
- Gaussian and multi‑Gaussian model definitions  
- Synthetic noisy dataset generation  
- Peak detection using `scipy.signal.find_peaks`  
- Automatic parameter initialization from detected peaks  
- Nonlinear curve fitting with `scipy.optimize.curve_fit`  
- BIC calculation for model comparison  
- Plotting of noisy data and fitted curves  

## Limitations
- when peaks are overlapping to a significant degree, peak estimates will veer off from the true peak
- noise estimates are hard coded to a wave from 0 < x < 20, so any data that varies from that window requires a manual change to the code
- currently only tests BIC against reduced number of peaks, never tests adding peaks to what is found
- only tests the peaks that are closest together

## To-Do
- add area and error calculations
