import numpy as np
import matplotlib.pyplot as plt

# units are SI system
PI = np.pi

# Independent parameters

l = 0.5                # half length of antenna/ quarter wavelength
c = 2.9979e8           # speed of light
mu = 4*PI*1e-7         # permeability of free space
N = 100                # number of sections in each antenna half
Im = 1                 # amplitude of current 
a = 0.01               # radius of antenna wire

# Dependent parameters

lamda = l*4.0          # wavelength of radiation
f = c/lamda            # frequency of radiation

k = (2*PI)/lamda       # wavenumber 
dz = l/(N)             # length of each current sample

TotalCurrentElements = 2*N+1
TotalPotentialElements = TotalCurrentElements


# Problem Statement 1

z  = np.linspace(-N,N,TotalCurrentElements)      # locations of all 2N+1 current elements
z = z*l/N

I = np.zeros(TotalCurrentElements)               # current vector corresponding to z locations
I[0] = 0
I[N] = Im
I[2*N] = 0

u = np.delete(z,[0,N,2*N])                       # locations of 2N-2 unknown current elements                   

J = np.zeros(2*N-2)                              # current vector corresonding to u locations


# Problem Statement 2

def matrixM(N):                                  # function to calculate M, where (2pi a)H = M*J
    M = np.zeros((2*N-2,2*N-2), dtype=float)
    np.fill_diagonal(M,1)
    return (1/(2*PI*a))*M

M = matrixM(N)


# Problem Statement 3

R = np.zeros((2*N+1,2*N+1),dtype=float)          # declaring a matrix of size 2N+1 by 2N+1 filled with zeros

'''
To avoid using nested for loop I used meshgrid method to get two matrices zy and zx 
such that zy varies horizontally from -0.5 to 0.5 while remains constant vertically,
similarly zx varies vertically from -0.5 to 0.5 while remains constant horizontally.
The difference of these two can directly be used to calculate R matrix.
'''

zy,zx = np.meshgrid(z,z)                         
R = np.sqrt(a**2 + (zx-zy)**2)

X = np.copy(R)                                   # this matrix is used to calculate Pb and P

X = np.delete(X,[0,N,2*N],axis=1)

R1n = X[N,:]                                     # extracting the R1n row to be used to
Pb = (((mu*dz)/(4*PI))*(np.exp(-1j*k*R1n)))/R1n  # calculate the value of Pb
                                                 # important to observe this done after deleting 
X = np.delete(X,[0,N,2*N],axis=0)                # columns but before deleting rows


P = np.zeros((2*N-2,2*N-2), dtype= complex)
P = ((mu*dz/(4*PI))*np.exp(-1j*k*X))/X


# Problem Statement 4

Q = np.zeros((2*N-2,2*N-2), dtype=complex)

Q = (a/mu)*P*(((1j*k)/X)+(1/(X**2)))             # calculating the Q and Qb matrices using the relations
                                                 # derived in the given assignment 
Qb = (a/mu)*Pb*((1j*k/R1n)+(1/(R1n**2)))        


# Problem Statement 5

J = np.dot(np.linalg.inv(M-Q),Qb)

I[1:N] = abs(J[0:N-1])                           # to obtain I from J we assign the first half of J to
I[N+1:-1] = abs(J[N-1:])                         # I index from 1 to N and 2nd half of J to I index N+1 to -1

I1 = np.sin(k*(l-z))                             # current using given formula

h = np.linspace(0,2*N,2*N+1)                     # creating the x axis for plot
h = (h-N)*l/N

plt.plot(0)                                                # plotting the graph using matplotlib library
plt.plot(h, I, 'b-', label='Current estimated')
plt.plot(h, I1, 'g-', label='Current using equation')
plt.xlabel("Length of Antenna (feeder at origin) (metres)")
plt.ylabel("Current Value (Amperes)")
plt.legend()
plt.grid(True)
plt.show()