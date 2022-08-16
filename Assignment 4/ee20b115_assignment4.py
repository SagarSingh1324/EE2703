from pylab import *                        # importing necessary libraries
from scipy import special as sp
from scipy import integrate as inte

# problem statement 1

X = linspace(-2*pi,4*pi,401)
X = X[:-1]

def fexp(t):                               # defining function that return exp(X) when X is vector/scalar
        return exp(t)

def fcoscos(t):                            # defining function that return cos(cos(X)) when X is vector/scalar
        return cos(cos(t))

Y = fexp(X%(2*pi))                         # expected function using fourier coefficients
Z = fcoscos(X%(2*pi))                      # expected function using fourier coefficients

Y2 = exp(X)
Z2 = cos(cos(X))

figure(1)
title("true plot and expected plot of y = exp(x)")
semilogy(X,Y, 'r.--')                                   # expected plot using fourier expansion
semilogy(X,Y2, 'b')                                     # true plot acccording to the function
xlabel("x-axis")
ylabel("y-axis")
legend(["expected plot","true plot"],loc='upper right')
grid(True)
show()

figure(2)
title("true plot and expected plot of y = cos(cos(x))")
plot(X,Z, 'r.--')                                       # expected plot using fourier expansion
plot(X,Z2, 'b')                                         # true plot acccording to the function
xlabel("x-axis")
ylabel("y-axis")
legend(["expected plot","true plot"],loc='upper right')
grid(True)
show()

# problem statement 2

a01 = inte.quad(lambda x: exp(x), 0, 2*pi)       # finding he a0 coefficient of exp(x) fourier expansion
a01 = a01[0]/(2*pi)

A1 = []
B1 = []

for n in range(1,26):                            # finding the 'a' coefficients of y = exp(x) and storing in A1 list
        func = lambda x: exp(x)*cos(n*x)
        intg = inte.quad(func, 0, 2*pi)
        A1.append(intg[0]/pi)

for n in range(1,26):                            # finding the 'b' coefficients of y = exp(x) and storing in B1 list
        func = lambda x: exp(x)*sin(n*x)
        intg = inte.quad(func, 0, 2*pi)
        B1.append(intg[0]/pi)

a02 = inte.quad(lambda x: cos(cos(x)), 0, 2*pi)   # finding he a0 coefficient of cos(cos(x)) fourier expansion
a02 = a02[0]/(2*pi)

A2 = []
B2 = []

for n in range(1,26):                             # finding the 'a' coefficients of y = cos(cos(x)) and storing in A2 list
        func = lambda x: cos(cos(x))*cos(n*x)
        intg = inte.quad(func, 0, 2*pi)
        A2.append(intg[0]/pi)

for n in range(1,26):                             # finding the 'b' coefficients of y = cos(cos(x)) and storing in B2 list
        func = lambda x: cos(cos(x))*sin(n*x)
        intg = inte.quad(func, 0, 2*pi)
        B2.append(intg[0]/pi)     

F1 = []                                      # making the vector (a0,a1,b1,....,a25,b25) with 
F1.append(a01)                               # y = exp(x)'s coefficient and storing as F1
m=0
n=0
for i in range (1,51):
        if(i%2!=0):
                F1.append(A1[m])
                m = m+1
        elif(i%2==0):
                F1.append(B1[n]) 
                n = n+1      

F2 = []                                     # making the vector (a0,a1,b1,....,a25,b25) with 
F2.append(a02)                              # y = cos(cos(x))'s coefficient and storing as F2
m=0
n=0
for i in range (1,51):
        if(i%2!=0):
                F2.append(A2[m])
                m = m+1
        elif(i%2==0):
                F2.append(B2[n]) 
                n = n+1      

# problem statement 3

La = []                             # creating a vector (a0,a1,.....,a24,a25) and storing as La
for i in range (0,26):
        La.append('a'+str(i))

Lb = []                             # creating a vector (b1,.....,b24,b25) and storing as Lb
for i in range (1,26):
        Lb.append('b'+str(i))   

k = 1
l = 0
L = []                              # creating a list (a0,a1,b1,...,a25,b25) by joining La and Lb and storing as L
L.append(La[0])
for i in range (1,51):
        if (i%2!=0):
                L.append(La[k])
                k = k+1
        else:
                L.append(Lb[l])
                l = l+1        

figure(3)                          # plotting the magnitude of fourier coefficients of exp(x) vs L vector using semilogy
semilogy(L,absolute(F1), 'r.')
title("magnitude of fourier coefficients of y = exp(x)")
ylabel("logarithmic y-axis")
xlabel("linear x-axis")
legend(["fourier coefficients"],loc='upper right')
grid(True)
show()
 
figure(4)                          # plotting the magnitude of fourier coefficients of exp(x) vs L vector  using loglog
loglog(L,absolute(F1), 'r.')
title("magnitude of fourier coefficients of y = exp(x)")
ylabel("logarithmic y-axis")
xlabel("logarithmic x-axis")
legend(["fourier coefficients"],loc='upper right')
grid(True)
show()

figure(5)                          # plotting the magnitude of fourier coefficients of cos(cos(x)) vs L vector using semilogy
semilogy(L,absolute(F2), 'r.')
title("magnitude of fourier coefficients of y = cos(cos(x))")
ylabel("logarithmic y-axis")
xlabel("linear x-axis")
legend(["fourier coefficients"],loc='upper right')
grid(True)
show()

figure(6)                          # plotting the magnitude of fourier coefficients of cos(cos(x)) vs L vector using loglog
loglog(L,absolute(F2), 'r.')
title("magnitude of fourier coefficients of y = cos(cos(x))")
ylabel("logarithmic y-axis")
xlabel("logarithmic x-axis")
legend(["fourier coefficients"],loc='upper right')
grid(True)
show()

# problem statement 4

L1 = linspace(0,2*pi,401)
L1 = L1[:-1]                          # drop last term to have a proper periodic integral

matA = zeros((400,51))                # creating a (400,51) matrix with all values initialized to 0
matA[:,0]=1 
for k in range(1,26):                 # filling up the matrix using L1 vector
        matA[:,2*k-1]=cos(k*L1) 
        matA[:,2*k]=sin(k*L1) 

# problem statement 5

b1 = exp(L1)                               # true values at L1 sample points for exp(x)
b2 = cos(cos(L1))                          # true values at L1 sample points for cos(cos(x))

C1 = lstsq(matA,b1, rcond=None)[0]         # the ’[0]’ is to pull out the best fit vector
C2 = lstsq(matA,b2, rcond=None)[0]         # the ’[0]’ is to pull out the best fit vector

figure(7)                                  # absolute values of coefficients using fourier and least sq. 
semilogy(L,absolute(F1), 'r.')           # using semilogy for function y = exp(x)
semilogy(L,absolute(C1), 'g+')
title("values of coeffs. using fourier vs using least square")
legend(["fourier coeffs.","least square coeffs."],loc='upper right')
ylabel("logarithmic y-axis")
xlabel("linear x-axis")
grid(True)
show()

figure(8)                                  # absolute values of coefficients using fourier and least sq. 
semilogy(L,absolute(F2), 'r.')           # using semilogy for function y = cos(cos(x))
semilogy(L,absolute(C2), 'g+')
title("values of coeffs. using fourier vs using least square")
legend(["fourier coeffs.","least square coeffs."],loc='upper right')
ylabel("logarithmic y-axis")
xlabel("linear x-axis")
grid(True)
show()

figure(9)                                  # absolute values of coefficients using fourier and least sq. 
loglog(L,absolute(F1), 'r.')               # using loglog for function y = exp(x)
loglog(L,absolute(C1), 'g+')
title("values of coeffs. using fourier vs using least square")
legend(["fourier coeffs.","least square coeffs."],loc='upper right')
ylabel("logarithmic y-axis")
xlabel("logarithmic x-axis")
grid(True)
show()

figure(10)                                 # absolute values of coefficients using fourier and least sq. 
loglog(L,absolute(F2), 'r.')               # using loglog for function y = cos(cos(x))
loglog(L,absolute(C2), 'g+')
title("values of coeffs. using fourier vs using least square")
legend(["fourier coeffs.","least square coeffs."],loc='upper right')
ylabel("logarithmic y-axis")
xlabel("logarithmic x-axis")
grid(True)
show()

# problem statement 6

C1F1 =[]                                     # defining C1F1 as difference of corresponding coefficients using least sq. and fourier 
for i in range(0,51):
        C1F1.append(C1[i]-F1[i])        

maxdev1 = max(absolute(C1F1))                # finding maximum deviation between least sq. and fourier for exp(x)
print("\nMaximun deviation of coefficient using least square vs using fourier is: "+str(maxdev1))

C2F2 =[]                                     # defining C1F1 as difference of corresponding coefficients using least sq. and fourier
for i in range(0,51):
        C2F2.append(C2[i]-F2[i])

maxdev2 = max(absolute(C2F2))                # finding maximum deviation between least sq. and fourier for cos(cos(x))
print("\nMaximun deviation of coefficient using least square vs using fourier is: "+str(maxdev2))

# problem statement 7
 
def func(t):                                  # defining the function that computes matrix A used  
        matA = zeros((400,51))                # in lstsq function for input vector 
        matA[:,0]=1 
        for k in range(1,26):
                matA[:,2*k-1]=cos(k*t) 
                matA[:,2*k]=sin(k*t) 
        return matA          

A11 = func(X)                                 # A11 is the x-axis from -2pi to 4 pi divided in 400 parts

prod1 = dot(A11,F1)                           # computing exp(x) values using coefficients estimated using fourier formula
prod2 = dot(A11,C1)                           # computing exp(x) values using coefficients estimated using lstsq method
prod3 = dot(A11,F2)                           # computing cos(cos(x)) values using coefficients estimated using fourier formula
prod4 = dot(A11,C2)                           # computing cos(cos(x)) values using coefficients estimated using lstsq  method

figure(11)
title("y = exp(x) using fourier and least square")
plot(X,prod1, 'r.--')                                      # plotting y = exp(x) using fourier vs x axis
plot(X,prod2, 'b')                                         # plotting y = exp(x) using least square vs x axis
xlabel("x-axis")                                           # formatting 
ylabel("y-axis")
legend(["fourier method","least square method"],loc='upper right')
grid(True)
show()      

figure(12)
title("y = cos(cos(x)) using fourier and least square")
plot(X,prod3, 'r.--')                                      # plotting y = cos(cos(x)) using fourier vs x axis
plot(X,prod4, 'b')                                         # plotting y = cos(cos(x)) using least square vs x axis
xlabel("x-axis")                                           # formatting
ylabel("y-axis")
legend(["fourier method","least square method"],loc='upper right')
grid(True)
show()

#   ---x---   #