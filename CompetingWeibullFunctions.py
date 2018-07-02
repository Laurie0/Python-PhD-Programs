# The system we are modelling consists of samples from which a crystal with 1 of 2 different
# structures (alpha or gamma) can form after some time t.

# This program fits two Weibull functions to experimental data in the form of the fractions 
# of samples which haven't yet crystallised in the form of the given crystal structure against time.

# Weibull functions have two parameters: beta and tau and take the form - P(t) = exp(- (t/tau)^Beta ) 

import numpy as np
from scipy.optimize import fmin

# ----------------------------------- hazard function---------------------------------------

# This generates the hazard function value of a Weibull distribution of a given beta and tau 
# at at a specified time, t

def hazard(t, beta, tau):

    if t == 0:
        output = 0
    else:
        output = beta * (t**(beta-1)) / (tau**beta)
    return output

# --------------------- Euler method used to plot alpha and gamma CIFs ---------------------

# Here we use the euler method to plot alpha and gamma CIFs of given Beta and Tau values and
# compare them to our experimental data.

def euler(parameters):

    tau_a, beta_a, tau_g, beta_g = parameters

    P = 1.0
    I1 = 0.0
    I2 = 0.0

    t = 0.0
    t_step = 0.01

    x = []
    y = []
    y1 = []
    y2 = []
    x_save = []
    y_save = []
    y1_save = []
    y2_save = []

    # Reading in experimental data from text file

    data = np.genfromtxt("250mgml_All_runs.txt")
    data = zip(*data)

    # Here we are rounding the time intervals to same number of dp as the time step of our
    # euler method

    for a in data[0]:
        b = str(round(a,2))
        x_save.append(b)

    while t < float(data[0][-1]):

        # Calculating the value of the alpha and gamma CIFs (I1 and I2) at t + t_step

        P = P -  (P * (hazard(t, beta_a, tau_a) + hazard(t, beta_g, tau_g))   ) * t_step
        I1 = I1 + P * hazard(t, beta_a, tau_a) * t_step
        I2 = I2 + P * hazard(t, beta_g, tau_g) * t_step
        t = t + t_step

        x.append(t)
        y.append(P)
        y1.append(I1)
        y2.append(I2)

        # Here we save the values of the calculated CIFs at the times at which our 
        # experimental data are recorded.

        for a in x_save:
            if str(t) == a:
                y_save.append(P)
                y1_save.append(I1)
                y2_save.append(I2)            

    # Here we calculate the square sum of the differences between our experimental data
    # and our generated functions.

    sum_sq = 0
    count = 0

    while count < len(data[0]):
        sum_sq = sum_sq + ((y1_save[count] - 1.0+data[1][count])**2)
        sum_sq = sum_sq + ((y2_save[count] - 1.0+data[2][count])**2)
        count += 1

    print parameters
    print sum_sq
    return sum_sq

# ------------------------- function for plotting with final params ------------------------

# This is almost identical to the previous function but here we plot the data and output the
# results to a text file.

def euler_plotter(parameters):

    tau_a, beta_a, tau_g, beta_g = parameters
    P = 1.0
    I1 = 0.0
    I2 = 0.0
    t = 0.0
    t_step = 0.01

    x = []
    y = []
    y1 = []
    y2 = []
    x_save = []
    y_save = []
    y1_save = []
    y2_save = []

    data_list = np.genfromtxt("250mgml_All_runs.txt")
    data_list = zip(*data_list)
    data=np.asarray(data_list)


    for a in data[0]:
        b = str(round(a,2))
        x_save.append(b)

    while t < float(data[0][-1]):

        P = P -  (P * (hazard(t, beta_a, tau_a) + hazard(t, beta_g, tau_g))   ) * t_step
        I1 = I1 + P * hazard(t, beta_a, tau_a) * t_step
        I2 = I2 + P * hazard(t, beta_g, tau_g) * t_step

        t = t + t_step

        for a in x_save:
            if str(t) == a:
                y_save.append(P)
                y1_save.append(I1)
                y2_save.append(I2)            

        x.append(t)
        y.append(P)
        y1.append(I1)
        y2.append(I2)

    # The minimized function data is written to a text file and plotted.

    import matplotlib.pyplot as plt
    
    f = open('CompetingWeibullFunctionFits.txt', 'w')

    for i in range(0,len(data[0])):
        data[1][i]=1.0-data[1][i]
        data[2][i]=1.0-data[2][i]
        f.write(str(data[0][i]) + '   '  + str(y1_save[i]) + '   ' + str(y2_save[i]) + '   ' + str(y_save[i]) +'\n')


    plt.plot(data[0], data[1], 'o')
    plt.plot(data[0], data[2], 'o')
    plt.plot(data[0], y1_save, '.')
    plt.plot(data[0], y2_save, '.')
    plt.plot(x, y, '.')
    plt.show()

# -------------------------------------- Main program --------------------------------------
        
# We start with some intial guess of parameters

guess = np.array([250.0, 1.5, 600.0, 0.5])

# We then use this guess as the starting point for the nelder-mead minimization algorithm

res = fmin(euler, guess, xtol=0.01, ftol=0.01)
print res


euler_plotter(res)

