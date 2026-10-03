import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

#returns the y-value of a gaussian peak based on inputted values
def gaussian (x, mu, sigma, amp):
    if sigma == 0:
        sigma = 1 * 10**-6
    else:
        return  (amp  *(np.exp( -((x-mu)**2) / (2*(sigma**2)))))

#takes n curves and returns parameters for curve fit
def n_gaussian (x, *params):
     result = 0
     for i in range(0, len(params), 3):
        result += gaussian(x,            # x value
                           params[i],    # mu
                           params[i+1],  # sigma
                           params[i+2])  # amp
     return result

#runs a for loop collecting data from gaussian w start value, end value, total num of data points
x = np.linspace(0,20,500)
true_curve = [10,0.7,8,6,1.6,6]
y_true = n_gaussian(x, *true_curve)

#adds random noise to simulate an experimental dataset
noise = 0.4
y_noisy = y_true + np.random.normal(0,noise, size=x.shape)

#manipulating constraints of find peaks to limit effect of noise
sigma_min = 1 # can be variable
sample_spacing = x[1] - x[0]

#noise estimate
edge_region = y_noisy[(x < 3) | (x > 17)]
noise_estimate = np.std(edge_region)
height_threshold = np.mean(edge_region) + 4 * noise_estimate

#identifies peaks in the noisy data and prints them
peaks_idx, properties = find_peaks(y_noisy, distance = int(sigma_min/sample_spacing), prominence = noise_estimate*3.5, height = height_threshold)
x_peaks = x[peaks_idx]
y_peaks = properties['peak_heights']
print(f"x peaks: {x_peaks}")
print(f"y peaks: {y_peaks}")

#prints the x and y values along a curve along w the noisy data
plt.scatter(x, y_noisy, label = "noisy signal")

popt = []
if len(peaks_idx) == 1:
    center = x[peaks_idx[0]]
    height = properties['peak_heights'][0]
    offset = 1  # variable; fails when peaks are sharp and close together
    n = len(x)

    # tests for single peak fit and error
    single_p0 = [center, 1, height]
    single_bounds = ([-np.inf, 0, 0], [np.inf, np.inf, np.inf])
    single_popt, single_pcov = curve_fit(gaussian,x,y_noisy,p0=single_p0,bounds = single_bounds)
    single_fit = gaussian(x, *single_popt)
    single_rss = np.sum((y_noisy - single_fit)**2)
    single_aic = n * np.log(single_rss/n) + 2*3
    print(f"single popt: {single_popt}")
    print(f"single aic: {single_aic}")

    #and compares to two peak fit
    double_p0 = [center - offset, 1, height/2, center + offset, 1, height/2]
    double_bounds = ([-np.inf, 0, 0, -np.inf, 0, 0],      # lower bounds
                     [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]) # upper bounds
    double_popt, double_pcov = curve_fit(n_gaussian, x, y_noisy, p0 = double_p0, bounds = double_bounds)
    double_fit = n_gaussian(x, *double_popt)
    double_rss = np.sum((y_noisy - double_fit)**2)
    double_aic = n * np.log(double_rss/n) + 2*6
    print(f"double popt: {double_popt}")
    print(f"double aic: {double_aic}")

    #sets popt equal to whichever peak has lower aic
    if single_aic + 10 > double_aic:
         popt,pcov = double_popt, double_pcov
    else:
         popt,pcov = single_popt, single_pcov

elif len(peaks_idx) > 1:
    p0 = []
    lower_bounds = []
    upper_bounds = []
    for i in range(0, len(peaks_idx), 1):
        p0.append(x[peaks_idx[i]])
        p0.append(1)
        p0.append(properties['peak_heights'][i])
        lower_bounds += [-np.inf,0,0]
        upper_bounds += [np.inf,np.inf,np.inf]
    bounds = (lower_bounds, upper_bounds)
    popt, pcov = curve_fit(n_gaussian,x,y_noisy,p0=p0,bounds=bounds)


#least-squares regression; fits the given p0 curve to the noisy data set and prints popt, the parameters for the fitted curve
print(f"true curve: {true_curve}")
print(f"popt: {popt}")
est_curve = n_gaussian(x, *popt)
plt.plot(x, est_curve, label = "fitted curve")
plt.xlabel("x")
plt.ylabel("approx distriibution")
plt.title("fitted curve")
plt.show()