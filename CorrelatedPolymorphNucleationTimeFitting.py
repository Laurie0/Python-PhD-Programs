# The system we are modelling consists of samples from which a crystal with 1 of 2 different
# structures (alpha or gamma) can form after some time t.

# This program fits two correlated functions to our experimental alpha and gamma crystallisation time data. We start with a 
# system of many samples. We assume each sample has a gamma crystallisation time and an alpha crystallisation time. 
# Both sets of crystallisation times are Weibull distributed and those two Weibull
# distributions have different Beta and Tau  parameters. As well as being Weibull distributed, the alpha and gamma 
# crystallisation times for each sample are correlated such that a sample with a short alpha crystallisation time is
# more likely to have a short gamma crystallisation time and vice versa. For each sample the crystallisation time,
# is whichever is shortest of the alpha crystallisation time and the gamma crystallisation time for that sample. The
# structure of the crystal that forms is hence also decided by the shorter crystallisation time.

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import random as r
import datetime

# ------------------------weibull distributed random number function------------------------

# Function takes in a uniform random random number (uniform_random), and converts it to a
# random variable from a Weibull distribution characterized by parameters beta and tau.

def rand_weibull(tau, beta, uniform_random):

    output = tau*(   (   -np.log(1.0-uniform_random)   )**(1.0/beta)   )
    return output

#---------------------------------- Monte Carlo Function -----------------------------------

def monte(parameters, samples):

    tau_1, beta_1, tau_2, beta_2 = parameters
    
    # Here we create 2 distributions x and y
    # They are normally distributed but correlated according to the covariance matrix 'cov'

    mean = [0,0]
    cov = [[1, 0.95], [0.95, 1]]  # diagonal covariance
    
    x, y = np.random.multivariate_normal(mean, cov, samples).T
    
    # Here we change the normally distributed numbers to uniformly distributed numbers 
    # between 0 and 1    

    x = ss.norm.cdf(x)
    y = ss.norm.cdf(y)
    
    # Here we change the uniformly distributed numbers to Weibull distributed numbers

    count = 0
    for a in x:
        x[count] = rand_weibull(tau_1, beta_1, a)
        count += 1

    count = 0
    for a in y:
        y[count] = rand_weibull(tau_2, beta_2, a)
        count += 1

    # We assume each sample has some characteristic alpha induction time and gamma
    # induction time. For each pair of induction times (alphatime and gammatime)
    # the shortest induction time is selected as the induction time of the sample.

    count = 0
    alpha_nucleation_times = []
    gamma_nucleation_times = []
    
    while count < samples:
    
        alpha_time = x[count]
        gamma_time = y[count]
    
    
        if alpha_time < gamma_time:
            alpha_nucleation_times.append(alpha_time)
        elif alpha_time > gamma_time:   
            gamma_nucleation_times.append(gamma_time)
        elif alphatime == gammatime:
            print 'equality', alpha_time, gamma_time
        count = count + 1   
    
    # Retrieving actuall experimental data from text file

    data = np.genfromtxt("250mgml_All_runs.txt")
    data = zip(*data)

    # Here, we convert the lists of induction times to the fraction of samples where an
    # alpha crystal hasn't yet nucleated as a function of time (P_alpha), and the fraction of
    # samples where a gamma crystal hasn't yet nucleated as a function of time (P_gamma).
    
    P_alpha = []
    for a in data[0]:
        count = 0.0
        for b in alpha_nucleation_times:
            if b < a:
                count += 1.0
        P_alpha.append(1.0-count/float(samples))
    
    P_gamma = []
    for a in data[0]:
        count = 0.0
        for b in gamma_nucleation_times:
            if b < a:
                count += 1.0
        P_gamma.append(1.0-count/float(samples))

    # Comparing Simulated P(t)'s to experimental data and finding the square sum difference
    # of the two functions.
    
    sum_sq = 0
    count = 0
    while count < 190: 
        sum_sq = sum_sq + ((P_gamma[count] - data[1][count])**2)
        sum_sq = sum_sq + ((P_alpha[count] - data[2][count])**2)
        count += 1

    return sum_sq
    
#----------------------------- Simulated annealing function --------------------------------

def annealing(parameters, samples, temperature, variance, iteration_threshold):

    prev_guess = 1000000.0
    best_guess = 1000000.0
    best_guess_parameters = 0
    count = 0.0
    tau1, b1, tau2, b2 = parameters

    # In this function we try to find the set of beta and tau parameters for the weibull 
    # distributions, that minimizes the sum square difference between our simulated function
    # and experimental data

    while True:
        
        # If a new minima has not been found for some number iterations 

        if count > iteration_threshold:
            print 'done'
            print 'best parameters', best_guess_parameters, 'best guess sum_sq', best_guess
            return best_guess_parameters
            break

        count += 1.0

        # Saving parameters from previous iteration of loop

        guess_save = parameters
    
        # Testing new parameters which are slightly different from thos tested in the  
        # previous iteration.
    
        tau1 = r.normalvariate(1, variance) * tau1
        tau2 = r.normalvariate(1, variance) * tau2
        b1 = r.normalvariate(1, variance) * b1
        b2 = r.normalvariate(1, variance) * b2
        parameters = np.array([tau1, b1, tau2, b2])

        # Putting new parameters into our monte carlo function

        new_guess = monte(parameters, samples)

        # Here we test if the new parameters give a better outcome than those in the 
        # previous iteration.

        if new_guess < prev_guess:
            
            # If the new parameters give a better result than the previous guess we use the
            # new parameters as the start point for the next iteration.

            prev_guess = new_guess

            # Here if the guess is better than the best guess so far, we also re-asign that.
            
            if new_guess < best_guess:
                best_guess = new_guess
                best_guess_parameters = parameters
                count = 0
        
        # If the new guess is worse than the previous guess we either go back to the
        # previous parameters or stick with the new parameters. This depends on whether a 
        # uniform random number is greater than the exponential function shown below
        # which looks at the difference of the previous and new sum square differences.

        elif r.uniform(0,1) > np.exp((1/temperature)*(prev_guess - new_guess)):
            tau1, b1, tau2, b2 = guess_save
    
        else:
            prev_guess = new_guess

# --------------------------------------- Main Program -------------------------------------

# Initial guess parameters

tau1 = 250.0
b1 = 0.88
tau2 = 150.0
b2 = 0.59
guess = np.array([tau1, b1, tau2, b2])

# Parameters for simulated annealing and monte carlo functions

variance = 0.03
numberofsamples = 10**2
iteration_threshold = 10**2

# We call the simulated annealing function, reducing the 'temperature' on each iteration

counter = 0.0
print datetime.datetime.now().time()

while counter <= 0.7:

    temperature = (0.8 - counter)* 2 
    guess = annealing(guess, numberofsamples, temperature, variance, iteration_threshold)
    print counter
    counter += 0.1
    print datetime.datetime.now().time()

# We use the simulated annealing one last time with a large sample size for the monte carlo
# function, and a large iteration threshold to get more precise minima parameter values.

variance = 0.01
numberofsamples = 10**4
iteration_threshold = 10**3
temperature = 0.2

guess = annealing(guess, numberofsamples, temperature, variance, iteration_threshold)

print datetime.datetime.now().time().datetime.now().time()

