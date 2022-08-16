from pylab import *
from matplotlib import *
from scipy.linalg import lstsq

figNum=0

def spectrum(figTitle, w, Y, magStyle='g.-', phaseStyle='r.'):
    suptitle(figTitle)
    subplot(211)
    grid()
    plot(w, abs(Y), magStyle, lw=2)
    subplot(212)
    grid()
    plot(w, angle(Y), phaseStyle, lw=2)
    show()

def signal(t, x, figTitle, style='g.-'):
    title(figTitle)
    grid()
    plot(t, x, style)
    show()

def hammingWindow(n):
    N = n.size
    window = zeros(N)
    window = 0.54 + 0.46*cos((2*pi*n)/(N-1))
    return fftshift(window)

# Example 1 - sin(sqrt(2)t) [Without windowing]
t = linspace(-pi, pi, 65)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
y = sin(sqrt(2)*t)
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
spectrum(r"Spectrum of $sin(\sqrt{2}t)$", w, Y)

#[Windowing with Hamming Window]
t = linspace(-pi, pi, 65)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
n = arange(64)
y = sin(sqrt(2)*t) * hammingWindow(n)
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65)[:-1]
spectrum(r"Spectrum of $sin(\sqrt{2}t) * w(t)$", w, Y)


# Problem Statement 2 - spectrum of (cos(0.86 t))**3 [Without windowing]
t = linspace(-4*pi, 4*pi, 257)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
y = cos(0.86*t)**3
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/256.0
w = linspace(-pi*fmax, pi*fmax, 257)[:-1]
spectrum(r"Spectrum of $cos^3(0.86t)$", w, Y)

#[Windowing with Hamming Window]
t = linspace(-4*pi, 4*pi, 257)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
n = arange(256)
y = (cos(0.86*t))**3 * hammingWindow(n)
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y))/256.0
w = linspace(-pi*fmax, pi*fmax, 257)[:-1]
spectrum(r"Spectrum of $cos^3(0.86t) * w(t)$", w,Y )


def estimateWandD(w, wo, Y, do, pow=2):
    wEstimate = sum(abs(Y)**pow * abs(w))/sum(abs(Y)**pow) 
    print("wo = {:.03f}\t\two Estimated value = {:.03f}".format(wo, wEstimate))
    t = linspace(-pi, pi, 129)[:-1]
    y = cos(wo*t + do)

    c1 = cos(wEstimate*t)
    c2 = sin(wEstimate*t)
    A = c_[c1, c2]
    vals = lstsq(A, y)[0]
    dEstimate = arctan2(-vals[1], vals[0])
    print("do = {:.03f}\t\tdo Estimated value = {:.03f}".format(do, dEstimate))


# Problem Statement 3 - Estimation of w, d in cos(wt + d)
wo = 1.35
d = pi/2

print("Question 3:")
t = linspace(-pi, pi, 129)[:-1]
trueCos = cos(wo*t + d)
fmax = 1.0/(t[1]-t[0])
n = arange(128)
y = trueCos.copy()*hammingWindow(n)
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
spectrum(r"Spectrum of $cos(\omega_o t + \delta) \cdot w(t)$", w, Y)
estimateWandD(w, wo, Y, d, pow=1.75)

# Problem Statment 4 - Estimation of w, d in noisy cos(wt + d)
print("\nQuestion 4:")
trueCos = cos(wo*t + d)
noise = 0.1*randn(128)
n = arange(128)
y = (trueCos + noise)*hammingWindow(n)
fmax = 1.0/(t[1]-t[0])
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax, pi*fmax, 129)[:-1]
spectrum(r"Spectrum of $(cos(\omega_o t + \delta) + noise) \cdot w(t)$", w, Y)
estimateWandD(w, wo, Y, d, pow=2.5)

# Problem Statement - DFT of chirped signal
def chirped(t):
    return cos(16*(1.5*t + (t**2)/(2*pi)))

t = linspace(-pi, pi, 1025)[:-1]
x = chirped(t)
signal(t, x, r"$cos(16(1.5 + \frac{t}{2\pi})t)$")
fmax = 1.0/(t[1]-t[0])
X = fftshift(fft(x))/1024.0
w = linspace(-pi*fmax, pi*fmax, 1025)[:-1]
spectrum(r"DFT of $cos(16(1.5 + \frac{t}{2\pi})t)$", w, X, 'g.-', 'r.')

n = arange(1024)
x = chirped(t)*hammingWindow(n)
signal(t, x, r" $cos(16(1.5 + \frac{t}{2\pi})t) \cdot w(t)$")
X = fftshift(fft(x))/1024.0
spectrum(r"DFT of $cos(16(1.5 + \frac{t}{2\pi})t) \cdot w(t)$", w, X, 'g.-', 'r.')

# Problem Statement 6 - Time evolution of DFT of chirped signal
def STFT(x, t, batchSize=64):
    t_batch = split(t, 1024//batchSize)
    x_batch = split(x, 1024//batchSize)
    X = zeros((1024//batchSize, batchSize), dtype=complex)
    for i in range(1024//batchSize):
        X[i] = fftshift(fft(x_batch[i]))/batchSize
    return X

# [plots of STFT]
def STFT3d(t, w, X, colorMap=cm.jet, showFig=True, blockFig=False):
    global figNum 
    t = t[::64]
    w = linspace(-fmax*pi,fmax*pi,65)[:-1]
    t, w = meshgrid(t, w)

    fig = plt.figure(figNum)
    ax = fig.add_subplot(211, projection='3d')
    surf = ax.plot_surface(w, t, abs(X).T, cmap=colorMap)
    fig.colorbar(surf)
    plt.xlabel(r"Frequency")
    plt.ylabel(r"Time")
    plt.title(r"Magnitude $\|Y\|$")

    ax = fig.add_subplot(212, projection='3d')
    surf = ax.plot_surface(w, t, angle(X).T, cmap=colorMap)
    fig.colorbar(surf)
    plt.xlabel(r"Frequency")
    plt.ylabel(r"Time")
    plt.title(r"Angle $\angle Y$")
    if showFig:
        plt.show(block=blockFig)
    figNum+=1

x = chirped(t)
X = STFT(x, t)
STFT3d(t, w, X)

x = chirped(t)*hammingWindow(arange(1024))
X = STFT(x, t)
STFT3d(t, w, X, colorMap=cm.brg, blockFig=True)

#########------xxxxxxxxxxxxxx------#########



