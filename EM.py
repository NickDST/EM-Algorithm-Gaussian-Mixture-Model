from sympy.interactive import printing
printing.init_printing(use_latex=True)
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


from sympy import Eq, solve_linear_system, Matrix
from numpy import linalg
import numpy as np
import sympy as sp
from sympy import *

# pi1, x, mu1, mu2, sigma1, sigma2, sigma,  w_n, u, v1, v2, theta, r, l, k = sp.symbols('pi1, x mu1 mu2 sigma1 sigma2 sigma w_n u v1, v2 theta r lambda, k')
x, sigma1, sigma2,  w_n, theta, r, l, k = sp.symbols('x sigma1 sigma2 w_n theta r lambda, k')

mu1 = 5
sigma = 5
mu2 = 30
pi1 = 0.4

num_observations = 100

### generate some normal variables from a given set of mu1, sigma1, mu2, sigma2

observations = []

for i in range(num_observations):
    latent_state = np.random.uniform(0,1,1)[0]
    if(latent_state > pi1):
        observations.append(np.random.normal(mu1, sigma))
    else:
        observations.append(np.random.normal(mu2, sigma))


pi1, mu1, mu2, sigma = sp.symbols('pi1 mu1 mu2 sigma')


gaussian1 = exp(-((x-mu1)**2 / (2 * sigma**2)))
gaussian2 = exp(-((x-mu2)**2 / (2 * sigma**2)))

gamma1 = (pi1 * gaussian1.subs(x, w_n)) / ((pi1 * gaussian1.subs(x, w_n) + (1 - pi1) * gaussian2.subs(x, w_n)))
gamma2 = ((1 - pi1) * gaussian2.subs(x, w_n)) / ((pi1 * gaussian1.subs(x, w_n) + (1 - pi1) * gaussian2.subs(x, w_n)))





### calculate pi1_new
def guess_pi1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess):

    numerator = 0
    denominator = 0
    for i in range(len(observations)):
        numerator += gamma1.subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)
        denominator += (gamma1+gamma2).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)

    new_pi1_guess = numerator / denominator
    return new_pi1_guess


### calculate u1_new
def guess_mu1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess):

    numerator = 0
    denominator = 0
    for i in range(len(observations)):
        numerator += (gamma1 * w_n).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)
        denominator += (gamma1).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)

    new_mu1_guess = numerator / denominator
    return new_mu1_guess


### calculate u2_new
def guess_mu2(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess):

    numerator = 0
    denominator = 0
    for i in range(len(observations)):
        numerator += (gamma2 * w_n).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)
        denominator += (gamma2).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)

    new_mu2_guess = numerator / denominator
    return new_mu2_guess


### calculate v1_new
def guess_sigma1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess):

    numerator = 0
    denominator = 0
    for i in range(len(observations)):
        numerator += (gamma1 * (w_n - mu1)**2).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)
        denominator += (gamma1).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)

    new_sigma_guess = numerator / denominator
    return new_sigma_guess


### calculate v1_new
def guess_sigma2(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess):

    numerator = 0
    denominator = 0
    for i in range(len(observations)):
        numerator += (gamma2 * (w_n - mu2)**2).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)
        denominator += (gamma2).subs(mu1, mu1_guess).subs(mu2, mu2_guess).subs(sigma, sigma_guess).subs(w_n, observations[i]).subs(pi1, pi1_guess)

    new_sigma_guess = numerator / denominator
    return new_sigma_guess



iterations = 10

pi1_guess = 0.7
mu1_guess = 30
mu2_guess = 6
sigma_guess= 9


for i in range(iterations):

    pi1_guess = guess_pi1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess)
    mu1_guess = guess_mu1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess)
    mu2_guess = guess_mu2(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess)
    variance = guess_sigma1(gamma1, gamma2, mu1_guess, mu2_guess, sigma_guess, observations, pi1_guess)
    sigma_guess = sqrt(variance)
    print("########## Curr Round: ", i)
    print("pi1_guess: ", pi1_guess)
    print("mu1_guess: ", mu1_guess)
    print("mu2_guess: ", mu2_guess)
    print("sigma_guess: ", sigma_guess)


print(">>>>>>> Final values: ")
print("pi1_guess: ", pi1_guess)
print("mu1_guess: ", mu1_guess)
print("mu2_guess: ", mu2_guess)
print("sigma_guess: ", sigma_guess)
