#!/usr/local/bin/python3

from scipy.signal import butter, lfilter, freqz

import numpy as np
import math
import matplotlib.pyplot as plt

ORDER=8
FC=150
SR=44100


# Chart 1: filter frequency response
#plt.figure(1)
#plt.clf()


#for ORDER in [2, 4, 8]:
#  b, a = butter(N=ORDER, Wn=FC, btype='low', analog=False, output='ba', fs=SR)
#  w, h = freqz(b=b, a=a, fs=SR, worN=512)
#  plt.plot(w, abs(h), label="order = %d" % ORDER)

#plt.xlabel('Frequency (Hz)')
#plt.ylabel('Gain')
#plt.grid(True)
#plt.legend(loc='best')
#plt.xlim(0, 2000)

# fetch float data from file
D = []
O = []
with open("inputs02.raw", "r") as inputs:
  for line in inputs.readlines():
    D += [float(line)]

with open("outputs02.raw", "r") as outputs:
  for line in outputs.readlines():
    O += [float(line)]

# Chart 2: signal and filtered one
T = 0.05
nsamples = 2205
t = np.arange(0, nsamples) / SR
plt.figure(2)
plt.clf()
plt.plot(t, D, label='Guitar signal')
#plt.plot(t, O, label='Filtered signal')

#for ORDER in [2, 4]:
#  b, a = butter(N=ORDER, Wn=FC, btype='low', analog=False, output='ba', fs=SR)
#  y = lfilter(b, a, D)
#  plt.plot(t, y, label='Filtered signal (Order %g )' % ORDER)

plt.xlabel('time (seconds)')
plt.grid(True)
plt.axis('tight')
plt.legend(loc='upper left')

# chart 3: ideal signal:
plt.figure(3)
plt.clf()
N = 82
x = np.arange(0,0.05,0.00002)   # start,stop,step
y = np.sin(2*math.pi*N*x)
for H in [2, 3, 4, 5, 6]:
  y += np.sin(H*SR/N*x)
plt.xlim(0, 0.05)
plt.plot(x,y)

b, a = butter(N=2, Wn=FC, btype='low', analog=False, output='ba', fs=SR)
yp = lfilter(b, a, y)
plt.plot(x, yp, label='Filtered signal (Order 2 )')

plt.show()