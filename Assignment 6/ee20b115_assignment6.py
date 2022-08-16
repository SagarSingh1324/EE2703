import numpy as np
import scipy.signal as sp
import matplotlib.pyplot as plt

# problem statemment 1
plt.figure(0)
plt.title("Plot of $\ddot{x} + 2.25x = e^{-0.5t}cos(1.5t)u(t)$")

t = np.linspace(0,50,500)

H1 = sp.lti([1],[1,0,2.25])

H2 = sp.lti([1,0.5],np.polymul([1,0,2.25],[1,1,2.5]))

t,x1 = sp.impulse(H2,None,t)

plt.plot(t,x1)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.xlim(0,50)
plt.grid(True)
plt.savefig("figure_0")
plt.show()

# problem statement 2
plt.figure(1)
plt.title("Plot of $\ddot{x} + 2.25x = e^{-0.05t}cos(1.5t)u(t)$")

H3 = sp.lti([1,0.05],np.polymul([1,0,2.25],[1,0.1,2.2525]))

t,x2 = sp.impulse(H3,None,t)

plt.plot(t,x2)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.xlim(0,50)
plt.grid(True)
plt.savefig("figure_1")
plt.show()

# problem statement 3
plt.figure(2)
plt.title("Plot of $\ddot{x} + 2.25x = e^{-0.05t}cos(wt)u(t)$")

t = np.linspace(0,50,500)
k = np.linspace(1.4,1.6,5)

for i in k:
    u = np.cos(i*t)*np.exp(-0.05*t)*np.heaviside(t,1)
    t,y,svec=sp.lsim(H1,u,t)
    plt.plot(t,y)
    
plt.grid(True)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.xlim(0,50)
plt.legend(["w=1.40","w=1.45","w=1.50","w=1.55","w=1.60"])
plt.savefig("figure_2")
plt.show()

# problem statement 4
plt.figure(4)
plt.title("Plot of x(t)/y(t) vs t")

t = np.linspace(0, 20, 200)

H4 = sp.lti([1, 0, 2], [1, 0, 3, 0])
t, x = sp.impulse(H4,None,t)

H5 = sp.lti([2], [1, 0, 3, 0])
t, y = sp.impulse(H5,None,t)

plt.plot(t, y)
plt.plot(t, x)

plt.xlabel("t") 
plt.ylabel("x(t)/y(t)")
plt.legend(["$x(t): \ddot{x}+(x-y)=0,\dot{x}(0)=0,{x}(0)=1$","$y(t): \ddot{y}+2(y-x)=0,\dot{y}(0)=0,{y}(0)=0$"],loc = "upper right")
plt.grid(True)
plt.xlim(0,20)
plt.savefig("figure_3")
plt.show()


# problem statement 5
plt.figure(5)
plt.title("Plot of |H(jw)| vs w")

L = 1e-6
C = 1e-6
R = 100

H6 = sp.lti([1], [L*C, R*C, 1])
w,S,phi = H6.bode()

plt.semilogx(w,S)
plt.xlabel("w")
plt.ylabel("|H(jw)| (in dB)")
plt.xlim(1e3,1e9)
plt.grid(True)
plt.savefig("figure_4")
plt.show()


# problem statement 6
plt.figure(6)
plt.title("Plot of \N{GREEK CAPITAL LETTER PHI}(H(jw)) vs w")

plt.semilogx(w,phi)
plt.xlabel("w")
plt.ylabel(" \N{GREEK CAPITAL LETTER PHI}(H(jw)) in \N{DEGREE SIGN}")
plt.xlim(1e3,1e9)
plt.grid(True)
plt.savefig("figure_5")
plt.show()

# problem statement 7(a)
plt.figure(7)
plt.title("Plot of $v_o(t)$ vs t")

w1  = int(1e3)
w2 = int(1e6)

t = np.linspace(0, 0.1, w2)

u = np.cos((w1)*t)*np.heaviside(t,1)-np.cos((w2)*t)*np.heaviside(t,1)

t,y,svec=sp.lsim(H6,u,t)

plt.plot(t,y)
plt.ylim(-1,1)
plt.xlim(0,0.00003)    
plt.xlabel("t")
plt.ylabel("$v_o(t)$")
plt.grid(True)
plt.savefig("figure_6")
plt.show()

# problem statement 7(b)
plt.figure(8)
plt.title("Plot of $v_o(t)$ vs t")

t2,y2,svec2=sp.lsim(H6,u,t)

plt.plot(t,y)
plt.xlim(0,0.01)  
plt.xlabel("t")
plt.ylabel("$v_o(t)$")
plt.grid(True)
plt.savefig("figure_7")
plt.show()