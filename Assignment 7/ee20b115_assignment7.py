import numpy as np
import sympy as sym
import scipy.signal as sig
import matplotlib.pyplot as plt

def SYMtoLTI(SYMFunc):
    num, den = SYMFunc.as_numer_denom()
    num = sym.Poly(num, s)
    den = sym.Poly(den, s)
    numCoeffs = num.all_coeffs()
    denCoeffs = den.all_coeffs()
    for i in range(len(numCoeffs)):
        x = float(numCoeffs[i])
        numCoeffs[i] = x
    for j in range(len(denCoeffs)):
        x = float(denCoeffs[j])
        denCoeffs[j] = x
    return numCoeffs, denCoeffs

# Problem Statement 1
R1 = 1e4
R2 = 1e4
C1 = 1e-9
C2 = 1e-9
G = 1.58
Vi = 1

s = sym.symbols('s')
A = sym.Matrix([[0, 0, 1, -1/G],
    [-1/(1+s*R2*C2), 1, 0, 0],[0, -G, G, 1],
    [-1/R1-1/R2-s*C1, 1/R2, 0, s*C1]])
b = sym.Matrix([0,0,0,-Vi/R1])
V = A.inv()*b

voNum, voDen = SYMtoLTI(V[3])
circuit1 = sig.lti(voNum, voDen)

w=np.linspace(1, 1e6, int(1e6))
w, mag, phase = sig.bode(circuit1, w)
fig1=plt.figure(1)
fig1.suptitle('Bode Plot of Transfer Function of Lowpass Filter')
plt.subplot(2,1,1)
plt.semilogx(w, mag)
plt.ylabel('$20log(\|H(j\omega)\|)$')
plt.subplot(2,1,2)
plt.semilogx(w, phase)
plt.xlabel(r'$\omega \ \to$')
plt.ylabel(r'$\angle H(j\omega)$')
plt.savefig("Figure_1")
plt.show()

# Step response
t = np.linspace(0, 0.1, int(1e6))
time, voStep = sig.step(circuit1, None, t)
plt.figure(2)
plt.plot(time, voStep)
plt.title('Step Response of Lowpass Filter')
plt.xlabel('$t\ \to$')
plt.ylabel('$V_o(t)\ \to$')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_2")
plt.show()

# Problem Statement 2
Vi = np.heaviside(t, 1)*(np.sin(2e3*np.pi*t)+np.cos(2e6*np.pi*t))

plt.figure(3)
plt.plot(t, Vi)
plt.title('$V_i(t)=(sin(2x10^3\pi t)+cos(2x10^6\pi t))u(t)$ to Lowpass Filter')
plt.xlabel('t')
plt.ylabel('Vi(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_3")
plt.show()

time, vOut, rest = sig.lsim(circuit1, Vi, t)
plt.figure(4)
plt.plot(time, vOut)
plt.title('$V_o(t)$ for $V_i(t)=(sin(2x10^3\pi t)+cos(2x10^6\pi t))u(t)$ for Lowpass Filter')
plt.xlabel('t')
plt.ylabel('Vo(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_4")
plt.show()

# Problem Statement 3
def HighPass(R1, R3, C1, C2, G, Vi):
    A = sym.Matrix([[0, -1, 0, 1/G],
                    [s*C2*R3/(s*C2*R3+1), 0, -1, 0],
                    [0, G, -G, 1],
                    [-s*C2-1/R1-s*C1, 0, s*C2, 1/R1]])

    b = sym.Matrix([0,
                    0,
                    0,
                    -Vi*s*C1])
    return (A.inv()*b)[3]

Vo = HighPass(10000, 10000, 1e-9, 1e-9, 1.58, 1)
voNum, voDen = SYMtoLTI(Vo)
circuit2 = sig.lti(voNum, voDen)

w, mag, phase = sig.bode(circuit2, w)
fig5=plt.figure(5)
fig5.suptitle('Bode Plot of Transfer Function of HighPass Filter')
plt.subplot(2,1,1)
plt.semilogx(w, mag)
plt.ylabel('$20log(\|H(j\omega)\|)$')
plt.subplot(2,1,2)
plt.semilogx(w, phase)
plt.xlabel(r'$\omega \ \to$')
plt.ylabel(r'$\angle H(j\omega)$')
plt.savefig("Figure_5")
plt.show()

vi = np.heaviside(t, 1)*(np.sin(2e3*np.pi*t)+np.cos(2e6*np.pi*t))

plt.figure(6)
plt.plot(t, vi)
plt.title('$V_i(t)=(sin(2x10^3\pi t)+cos(2x10^6\pi t))u(t)$ to HighPass Filter')
plt.xlabel('t')
plt.ylabel('Vi(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_6")
plt.show()

time, vOut, rest = sig.lsim(circuit2, vi, t)
plt.figure(7)
plt.plot(time, vOut)
plt.title('$V_o(t)$ for $V_i(t)=(sin(2x10^3\pi t)+cos(2x10^6\pi t))u(t)$ for HighPass Filter')
plt.xlabel('t')
plt.ylabel('Vo(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_7")
plt.show()

# Problem Statement 4
ViDampedHighFreq = np.heaviside(t, 1)*(np.sin(2*np.pi*t))*np.exp(-t)
ViDampedLowFreq = np.heaviside(t, 1)*(np.sin(2e5*np.pi*t))*np.exp(-t)

# Output for Low Frequency
time, VoutDampedLowFreq, rest = sig.lsim(circuit2, ViDampedHighFreq, t)
plt.figure(8)
plt.plot(time, VoutDampedLowFreq)
plt.title('$V_o(t)$ for $V_i(t)=sin(2\pi t)e^{-t}u(t)$ for HighPass Filter')
plt.xlabel('t')
plt.ylabel('Vo(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_8")
plt.show()

# Output for High Frequency
time, VoutDampedHighFreq, rest = sig.lsim(circuit2, ViDampedLowFreq, t)
plt.figure(9)
plt.plot(time, VoutDampedHighFreq)
plt.title('$V_o(t)$ for $V_i(t)=sin(2x10^5\pi t)e^{-t}u(t)$ for HighPass filter')
plt.xlabel('t')
plt.ylabel('Vo(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_9")
plt.show()

# Problem Statement 5
time, voStep = sig.step(circuit2, None, t)
plt.figure(10)
plt.plot(time, voStep)
plt.title('Step Response of HighPass Filter')
plt.xlabel('t')
plt.ylabel('Vo(t)')
plt.xlim(0, 1e-3)
plt.grid(True)
plt.savefig("Figure_10")
plt.show()