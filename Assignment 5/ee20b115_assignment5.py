import numpy
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import sys

if (len(sys.argv) == 6):                          # checking if correct no. of arguments are provided
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    radius = int(sys.argv[3])
    Niter = int(sys.argv[4])
    tol = float(sys.argv[5])
else:
    print("Insufficient arguments provided")  
    sys.exit(0)  

x = array(linspace(int(-Nx/2), int(Nx/2), Nx))
y = array(linspace(int(Ny/2), int(-Ny/2), Ny))

phi = zeros([Nx,Ny])                              # making the grid over which contour plots will be made
[Y,X] = meshgrid(x,y)

Z = where((X*X + Y*Y) <= (radius)**2)             # defing the area of 1 volt potential
phi[Z] = 1.0                                      # marking the potential of rod area
phi1 = phi.copy()

Err = zeros([Niter])                               
for k in range(Niter): 
    oldphi = phi.copy()
    phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])     # algorithm to perform iterations
    phi[0,1:-1] = phi[1,1:-1]                                                           # setting boundary conditions
    phi[1:-1,0] =  phi[1:-1,1]
    phi[1:-1,-1] = phi[1:-1,-2]
    phi[-1,1:-1] = 0
    phi[Z] = 1.0                                                                        # resetting the central region potential after an iteration
    Err[k] = (np.max(absolute(phi-oldphi)))                                             # recording the error

def func_fit(S,c):                                              # defining function to find the A and B values of y=Aexp(Bx) fit
    A = zeros([S.size, 2])                                      # and also find the function values for that fit
    A[:,0] = 1
    A[:,1] = S
    b = lstsq(A, log(c), rcond = None)[0]                       # taking the first fit
    fit = exp(dot(A,b))
    return b, fit

Err2 = []                                                       # new array to only plot every 50th element
for i in range(int(Niter/50)):
    Err2.append(Err[50*i]) 

# plot of error vs iterations using semilog
figure(0)
title("semilog plot of error after "+str(Niter)+" iterations")
L1 = linspace(1,len(Err2),len(Err2)) 
semilogy(50*L1,Err2,'b.-')
xlabel("iterations")
ylabel("error")
legend(["error"], loc = 'upper right')
grid(True)
savefig("figure_0")
show()

# plot of error vs iterations using loglog
figure(1)
title("loglog plot of error after "+str(Niter)+" iterations")
loglog(50*L1,Err2,'b.-')
xlabel("iterations")
ylabel("error")
legend(["error"], loc = 'upper right')
grid(True)
savefig("figure_1")
show()


Err3 = np.delete(Err2,[0])                             # deleteing the error of first 100 iterations 
L2 = np.delete(L1,[0])                                 # as the loglog plot is linear till first 100 iterations

# plot of true error and models made using lstsq
figure(2)
title("semilog plot of true, fit for all, fit for values excluding first 100")
semilogy(50*L1,Err2,'g')
b1, fit1 = func_fit(L1,Err2)
b2, fit2 = func_fit(L2,Err3)
semilogy(50*L1, fit1, 'b.--')
semilogy(50*L2, fit2, 'r.--')
legend(["true value of error", "fit for all error values", "fit after excluding first 100 error values"], loc = "upper right")
xlabel("iterations")
ylabel("error")
grid(True)
savefig("figure_2")
show()

x = linspace(1,Niter,Niter)
b3, fit3 = func_fit(x, Err)

def cumError(N, A, B):                                     # function to calculate cumulative error
    return -(A/B)*exp(B*(N+0.5))

def stopCon(Niter, tol):                                   # function to the no. of iterations at the specified tolerance level
    cumErr = []
    for n in range(1, Niter):
        cumErr.append(cumError(n, exp(b3[0]), b3[1]))
        if(absolute(cumErr[n-1]) <= tol):
            return cumErr[-1], n     
    return cumErr[-1], Niter

finErr, nStop = stopCon(Niter, tol)
print("The stopping condition for given tolerance is %g iterations. At stopping condition, the error is %g" %(nStop, finErr))

# surface plot of phi after iterations
figure(3) 
title("3D surface plot of $\phi$ after "+str(Niter)+" iterations")
ax = axes(projection ='3d')
ax.plot_surface(Y, X, phi, cmap ='jet')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')
savefig("figure_3")
show()

# contour plot of phi before iterations
figure(4)
title("contour Plot of $\phi$ before iterations")
contourf(Y, X, phi1, cmap = cm.jet)
grid(True)
savefig("figure_4")
show()

# contour plot of phi after iteratons
figure(5)
title("contour Plot of $\phi$ after "+str(Niter)+" iterations")
contourf(Y, X, phi, cmap = cm.jet)
grid(True)
savefig("figure_5")
show()

# plot of current vectors
figure(6)                                                          
title("plot of current vectors and 1V region")
Jy = (phi[2:,1:-1] - phi[0:-2,1:-1])
Jx = (phi[1:-1,0:-2] - phi[1:-1,2:])
quiver( Y[1:-1,1:-1], X[1:-1,1:-1], Jx[:,::-1], Jy[:,::-1])
scatter(Y, X, phi1, color = 'Red', linewidths=5)
legend(["current", "1V region"])
grid(True)
savefig("figure_6")
show()





