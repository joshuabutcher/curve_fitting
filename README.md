# Multi‑Peak Gaussian Fitting with BIC Model Selection

<img width="639" height="462" alt="curve_fit_ex" src="https://github.com/user-attachments/assets/d3e94371-c7d0-4c49-b719-c38beed218fd" />


This project performs multi‑peak Gaussian fitting on noisy data using SciPy’s nonlinear least‑squares optimizer. It detects peaks, builds initial parameter guesses, fits one or more Gaussian components, and selects from different peak models using BIC (Bayesian Information Criterion).

## Features
- Gaussian and multi‑Gaussian model definitions  
- Synthetic noisy dataset generation  
- Peak detection using `scipy.signal.find_peaks`  
- Automatic parameter initialization from detected peaks  
- Nonlinear curve fitting with `scipy.optimize.curve_fit`  
- BIC calculation for model comparison  
- Plotting of noisy data and fitted curves
- Calculation of peak area and error propagation

## Limitations
- when peaks are overlapping to a significant degree, peak estimates will veer off from the true peak
- noise estimates are hard coded to a wave from 0 < x < 20, so any data that varies from that window requires a manual change to the code
- currently only tests BIC against reduced number of peaks, never tests adding peaks to what is found
- only tests manipulating the peaks that are closest together, and uses the amplitude of the taller peak which can result in error
- hidden peaks not detected by find_peaks are never tested for

## To-Do
- test BIC peak, area, and error calculations over many different seeds
- add a calibration curve and error propagation
- test for peaks not identified by find_peaks
