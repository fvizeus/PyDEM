import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

phi = 0.4
alpha_mean = 0.1
alpha_std = 0.01

def discrete_normal(zmin, zmax, n):
    z = np.linspace(zmin, zmax, n)
    x = np.empty(n + 1)
    x[0] = np.exp(-z[0]**2/2.0)/norm.cdf(z[0])
    x[-1] = np.exp(-z[-1]**2/2.0)/(1.0 - norm.cdf(z[-1]))
    x[1:-1] = (np.exp(-z[:-1]**2/2.0) - np.exp(-z[1:]**2/2.0))/(norm.cdf(z[1:]) - norm.cdf(z[:-1]))
    x /= np.sqrt(2.0*np.pi)

    p = np.empty(len(z) + 1)
    p[0] = norm.cdf(z[0])
    p[-1] = 1.0 - norm.cdf(z[-1])
    p[1:-1] = norm.cdf(z[1:]) - norm.cdf(z[:-1])
    
    return x, p

alphas, phis = discrete_normal(-4, 4, 81)
alphas = alphas*alpha_std + alpha_mean
phis = phi*phis

print sum(phis)

plt.subplot(121)
plt.vlines(alphas, 0, phis)

plt.subplot(122)
alphas_pdf = np.linspace(-5, 5, 10001)*alpha_std + alpha_mean
plt.plot(alphas_pdf, norm.pdf(alphas_pdf, loc=alpha_mean, scale=alpha_std))

plt.show()
