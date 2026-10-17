import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

#returns the y-value of a gaussian peak based on inputted values
def gaussian(x, mu, sigma, amp):
    if sigma <= 1e-6:
            sigma = 1e-6
    return  (amp  *(np.exp( -((x-mu)**2) / (2*(sigma**2)))))

# takes the area under a single gaussian peak
def peak_area(sigma, amp):
    area = float(amp * sigma * np.sqrt(2*np.pi))
    da_dsig = float(amp * np.sqrt(2*np.pi))
    da_damp = float(sigma * np.sqrt(2*np.pi))
    return area, da_dsig, da_damp

# takes n curves and returns an array with the areas under each curve, derivative of area with respect to sigma, and derivative of area with respect of amp
def n_peak_area(popt):
    area_array = []
    dads_array = []
    dada_array = []
    for i in range(0,len(popt),3):
        area, da_dsig, da_damp = peak_area(popt[i+1],popt[i+2])
        area_array.append(area)
        dads_array.append(da_dsig)
        dada_array.append(da_damp)

    return area_array, dads_array, dada_array

#takes n curves and returns parameters for curve fit
def n_gaussian(x, *params):
     result = 0
     for i in range(0, len(params), 3):
        result += gaussian(x,            # x value
                           params[i],    # mu
                           params[i+1],  # sigma
                           params[i+2])  # amp
     return result

#runs a for loop collecting data from gaussian with start value, end value, total num of data points
x = np.linspace(0,20,500)
true_curve = [12,0.5,8,  7,1.5,6,  15,0.2,4]
y_true = n_gaussian(x, *true_curve)

#adds random noise to simulate an experimental dataset
noise = 0.4
np.random.seed(0)
y_noisy = y_true + np.random.normal(0,noise, size=x.shape)

n = len(x)

def fit_from(mus, amps):
    #takes a list of mus and amps and pairs them before unpacking with a flat sigma / width of 1
    p0 = []
    for mu, h in zip(mus, amps):
        for v in (mu, 1, h):
            p0.append(v)
    peak_count = len(mus)
    bounds = ([x.min(), 0, 0]*peak_count, [x.max(),np.inf,np.inf]*peak_count) # bounds from x-min to x-max for mu and 0 to inf for sigma and amp
    popt, pcov = curve_fit(n_gaussian, x, y_noisy, p0=p0, bounds=bounds)
    rss = np.sum((y_noisy - n_gaussian(x, *popt))**2)
    bic = n*np.log(rss/n) + 3*peak_count*np.log(n)
    return popt, pcov, bic


#manipulating constraints of find peaks to limit effect of noise
sigma_min = 1 # can be variable
sample_spacing = x[1] - x[0]

#noise estimate | works for this data, must find flat area for this test if real data doesn't fit to 0<x<20 as well as this
low_baseline = 3
high_baseline = 17
edge_region = y_noisy[(x < low_baseline) | (x > high_baseline)]
noise_estimate = np.std(edge_region)
height_threshold = np.mean(edge_region) + 4 * noise_estimate

prominence = noise_estimate * 1.5 # let find_peaks over-detect peaks so BIC test can run (changed from prominence = noise_estimate * 3.5)

#identifies peaks in the noisy data and prints them
peaks_idx, properties = find_peaks(y_noisy, distance = int(sigma_min/sample_spacing), prominence = prominence, height = height_threshold)
x_peaks = x[peaks_idx]
y_peaks = properties['peak_heights']
print(f"x peaks: {x_peaks}")
print(f"y peaks: {y_peaks}")

#prints the x and y values along a curve along w the noisy data
plt.scatter(x, y_noisy, label = "noisy signal")
plt.xlabel("x")
plt.ylabel("probability density")
plt.title("gaussian distribution ")
plt.legend()

if len(x_peaks) < 1:
    raise Exception("no peaks detected")

if len(peaks_idx) > 1:
    # generalizes BIC error tests to all multi-peak functions
    x_peaks = list(x_peaks)
    y_peaks = list(y_peaks)

    while True:
        popt, pcov, bic_2 = fit_from(x_peaks, y_peaks) # models based on peaks given by find_peaks
        if len(x_peaks) == 1:
            break
        j = int(np.argmin(np.diff(x_peaks)))
        x_red = x_peaks.copy()
        y_red = y_peaks.copy()
        if y_peaks[j] < y_peaks[j+1]:
            x_red.pop(j)
            y_red.pop(j)
        else:
            x_red.pop(j+1)
            y_red.pop(j+1)
        
        popt1, pcov1, bic_1 = fit_from(x_red, y_red)
        if bic_1 < bic_2:
            x_peaks, y_peaks = x_red, y_red
        else:
            break

if len(x_peaks) == 1:
    center, height, offset = x_peaks[0], y_peaks[0], 1
    single_popt, single_pcov, single_bic = fit_from([center], [height])
    double_popt, double_pcov, double_bic = fit_from([center - offset, center + offset], [height/2, height/2])

    #sets popt equal to whichever peak has lower bic
    if single_bic > double_bic:
         popt,pcov = double_popt, double_pcov
    else:
         popt,pcov = single_popt, single_pcov

# calculates error by pulling from the area and derivatives returned by n_peak_area and variances given by pcov
area, dads, dada = n_peak_area(popt)
error_array = []
for i in range(0, len(popt), 3):
    # retrieves the variances of sigma and amplitude and covariance between the two from pcov
    var_sig = pcov[i+1, i+1]
    var_amp = pcov[i+2, i+2]
    cov_sa = pcov[i+1, i+2]
    k = int(i/3)
    peak_variance = dads[k]**2 * var_sig + dada[k]**2 * var_amp + 2 * dads[k] * dada[k] * cov_sa
    peak_error = float(np.sqrt(peak_variance))
    error_array.append(peak_error)



#least-squares regression; fits the given p0 curve to the noisy data set and prints all the calculations
print(f"true curve: {true_curve}")
print(f"popt: {popt}")
print(f"peak areas: {area}")
print(f"errors: {error_array}")

est_curve = n_gaussian(x, *popt)
plt.plot(x, est_curve, label = "fitted curve")
plt.xlabel("x")
plt.ylabel("approx distriibution")
plt.title("fitted curve")
plt.show()