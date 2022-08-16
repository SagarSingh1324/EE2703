from pylab import *
from matplotlib import *

plt.figure(1)
plt.suptitle(r'FFT of $sin(5t)$')
x = linspace(0,2*pi,129)
x = x[:-1]

w = linspace(-64,64,129)
w=w[:-1]

y = sin(5*x)

Y = fftshift(fft(y))/128

subplot(2,1,1)
plot(w,abs(Y),'g.-')
xlim([-10,10])
grid(True)

subplot(2,1,2)
plot(w,angle(Y),'r.')

xlim([-10,10])
grid(True)
show()

plt.figure(2)
plt.suptitle(r'AM Modulation with $(1+0.1cos(t))cos(10t)$')
x = linspace(0,8*pi,513)
x = x[:-1]

w = linspace(-64,64,513)
w=w[:-1]

y = (1+0.1*cos(x))*cos(10*x)

Y = fftshift(fft(y))/512

subplot(2,1,1)
plot(w,abs(Y),'g.-')
xlim([-20,20])
grid(True)

subplot(2,1,2)
plot(w,angle(Y),'r.')

xlim([-20,20])
grid(True)
show()

plt.figure(3)
plt.suptitle(r'Spectrum of $sin^3(t)$')
x = linspace(0,2*pi,129)
x = x[:-1]

w = linspace(-64,64,129)
w=w[:-1]

y = (sin(x))**3

Y = fftshift(fft(y))/128

subplot(2,1,1)
plot(w,abs(Y),'g.-')
xlim([-10,10])
grid(True)

subplot(2,1,2)
plot(w,angle(Y),'r.')

xlim([-10,10])
grid(True)
show()

plt.figure(4)
plt.suptitle(r'Spectrum of $cos^3(t)$')
x = linspace(0,2*pi,129)
x = x[:-1]

w = linspace(-64,64,129)
w=w[:-1]

y = (cos(x))**3

Y = fftshift(fft(y))/128

subplot(2,1,1)
plot(w,abs(Y),'g.-')
xlim([-10,10])
grid(True)

subplot(2,1,2)
plot(w,angle(Y),'r.')

xlim([-10,10])
grid(True)
show()

plt.figure(5)
plt.suptitle(r'Spectrum of $cos(20t + 5cos(t))$')
x = linspace(0,2*pi,129)
x = x[:-1]

w = linspace(-64,64,129)
w=w[:-1]

y = cos(20*x + 5*cos(x))

Y = fftshift(fft(y))/128

subplot(2,1,1)
plot(w,abs(Y),'g.-')
xlim([-40,40])
grid(True)

i = where(abs(Y)>1e-3)
subplot(2,1,2)
plot(w[i],angle(Y[i]),'r.')
xlim([-40,40])
grid(True)
show()

plt.figure(6)
plt.suptitle(r'Comparison of spectrum of $e^{-\frac{t^2}{2}}$')
t =  np.linspace(-8*pi, 8*pi, 1025)
t = t[:-1]
xTrueGaussian = np.exp(-(t**2)/2)
Y = fftshift(fft(ifftshift(xTrueGaussian)))*8/1024.0

YMag = np.abs(Y)
YPhase = np.angle(Y)
absentFreqs = np.where(YMag < 1e-3)
YPhase[absentFreqs] = 0
w = np.linspace(-40, 40, 1025)
w = w[:-1]
subplot(221)
plot(w, YMag,'g.-')
xlim([-10, 10])
title("Estimated Spectrum")
grid()
subplot(223)
plot(w, YPhase, 'r.')
xlim([-10, 10])
grid()

trueY = np.exp(-(w**2)/2)/np.sqrt(2*pi)
trueYMag = np.abs(trueY)
trueYPhase = np.angle(trueY)
subplot(222)
plot(w, trueYMag,'g.-')
xlim([-10, 10])
title("True Spectrum")
grid()
subplot(224)
plot(w, trueYPhase, 'r.')
xlim([-10, 10])

grid()
show()

meanError = np.mean(np.abs(trueY - Y))
print(r'Magnitude of mean error between actual and computed values of the Gaussian: ', meanError) 