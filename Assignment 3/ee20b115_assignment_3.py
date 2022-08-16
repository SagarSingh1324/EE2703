from pylab import *                        # importing necessary libraries
from scipy import special as sp


# Downloading the script and running it to generate a
# set of data. The data is written out to a file whose name is fitting.dat.


# Loading fitting.dat. The data consists of 10 columns. The first column is
# time, while the remaining columns are data. Extracting the raw data columns and storing in curve_data.

data = np.loadtxt("fitting.dat")

curve_data = [[],[],[],[],[],[],[],[],[],[]]

for i in range(len(data)):
    for j in range(len(data[i])):
        curve_data[j].append(data[i][j])

del curve_data[0]      #  deleting first element(list) which has time points as they are not needed

# Plotting the curves in Figure 0 and adding labels to indicate the amount of noise(with stdev sigma) in each. 
t = linspace(0,10,101)                            # creating the time points
list2 = around(logspace(-1,-3,9),3)               # creating sigma values and rounding off to 3 decimals

figure(0)
for i in range(len(curve_data)):
    plot(t,curve_data[i])

def g(t, A, B):                                   # defining the function to get true value
    return A*sp.jn(2,t) + B*t
A = 1.05
B = -0.105
trueFunction = g(t, A, B)

plot(t, trueFunction)    

title("Q4: Data to be fitted to theory")
xlabel("t")
ylabel("f(t) + noise($\sigma_{n}$)")    
grid()

# labelling the curves using legend
legend(['$\sigma_{} = {}$'.format(0, list2[0]),'$\sigma_{} = {}$'.format(1, list2[1]),'$\sigma_{} = {}$'.format(2, list2[2]),
'$\sigma_{} = {}$'.format(3, list2[3]),'$\sigma_{} = {}$'.format(4, list2[4]),'$\sigma_{} = {}$'.format(5, list2[5]),
'$\sigma_{} = {}$'.format(6, list2[6]),'$\sigma_{} = {}$'.format(7, list2[7]),'$\sigma_{} = {}$'.format(8, list2[8]),'true value'])
savefig('figure_0')
show()      

# Generating a plot of the first column of data with error bars. Plotting every 5th data item to make the plot
# readable. Also plotting the exact curve to see how much the data diverges.
figure(1)
xlabel('time')
title('Q5: Data points for $\sigma$ = 0.10 along with exact function')
plot(t, trueFunction,color='#000000')

errorbar(t[::5], curve_data[0][::5], 0.1, fmt='co')                     # plotting the error bars with 5th data point + formatting
legend(['f(t)','error bar'])
grid()
savefig('figure_1')
show()

# Constrcting M and p matrices and the G as M*p and seeing it is equal to trueFunction
jColumn = sp.jn(2,t)
M = c_[jColumn, t]
p = array([A, B])
G = dot(M,p)

MatProduct = c_[t,G]                          # created using product of matrices
actual = c_[t, trueFunction]                  # created using trueFunction


# Calculating the “mean squared error” between the data (curve_data) and the assumed model(trueFunction). 
# Using the first column of curve_data for this.
A = linspace(0,2,21)
B = linspace(-0.2,0,21)

epsilon = empty((len(A), len(B)))
for i in range(len(A)):
    for j in range(len(B)):
        epsilon[i][j] = mean(square(curve_data[0][:] - g(t[:], A[i], B[j])))

# Plotting the contour plot of epsilon(i,j)
figure(2)
cont_plot=contour(A,B,epsilon, levels=21)
xlabel("A")
ylabel("B")
title("Q8: contour plot of $\epsilon_{ij}$")
clabel(cont_plot, inline=1, fontsize=5)                   # formatting 
plot([1.05], [-0.105], 'bo')                              # plotting the true value
annotate("exact location", (1.05, -0.105))                # annotating the true value
grid()
savefig('figure_2')
show()

# Using the Python function lstsq from scipy.linalg to obtain the best estimate of A and B. Using the
# array created in part 6 as the argument. 
p, *rest = lstsq(M,trueFunction,rcond=None)

# Plotting the error in the estimate of A and B for
# different data files versus the noise with stdev sigma as figure 3.
perr = empty((9, 2))                                                 # creating an empty array
for k in range(len(curve_data)):                                     # looping through each column
    perr[k], *rest = lstsq(M, curve_data[k], rcond=None)   
Aerr = array([square(x[0]-p[0]) for x in perr])
Berr = array([square(x[1]-p[1]) for x in perr])
plot(list2, Aerr, '--yo')                                            # plotting + formatting
plot(list2, Berr, '--mo')
xlabel("noise standard deviation")
title("Q10: variation of error with noise")
ylabel("MS error")
legend(["Aerr","Berr"])
grid()
savefig('figure_3')
show()

# Replotting the above curves using loglog as figure 4. 
figure(4)
loglog(list2, Aerr, 'yo')                                            # plotting + colour
loglog(list2, Berr, 'mo')
errorbar(list2, Aerr, std(Aerr), fmt='yo')
errorbar(list2, Berr, std(Berr), fmt='mo')
xlabel("$\sigma_{n}$")                                               # formatting
title("Q11: variation of error with noise")
ylabel("MS error")
legend(["Aerr","Berr"])
grid()
savefig('figure_4')
show()






