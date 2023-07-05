#!/usr/local/bin/python3

import wave
import struct

from scipy.fft import fft, fftfreq
from scipy.signal import butter, lfilter, freqz, sosfilt
import scipy.fftpack
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

ORDER=4
FC=82*2*1.25
FILE="guitar_E.wav"

mpl.style.use('dark_background');

fig1, axs1 = plt.subplots(2, 1, figsize=(7, 10))

# Print original signal
with wave.open(FILE, "rb") as wavfile:
  astr = wavfile.readframes(wavfile.getnframes())
  FR = wavfile.getframerate();
  N = wavfile.getnframes()
  signal = struct.unpack("%ih" % (wavfile.getnframes()* wavfile.getnchannels()), astr)
  signal = [float(val) / pow(2,15) for val in signal]

t = np.arange(0, N) / FR
axs1[0].plot(t, signal, label='Guitar signal', color='C3')
axs1[0].set_title("Raw guitar signal waveform 82 Hz")
axs1[0].grid(color='grey', linewidth=1, linestyle='--',)

print("read "+str(len(signal))+" frames")
print("in the range "+str(min(signal))+" to "+str(max(signal)))
print("Frame rate:"+ str(FR))
print("# frames :"+ str(N))


# Print dfft spectrum
yf = np.fft.fft(signal)
xf = np.fft.fftfreq(len(yf), 1/FR)[:N//2]

axs1[1].plot(xf, np.abs(yf[0:N//2]), color='C5')
axs1[1].axis([0, 2000, 0, 17000])
axs1[1].set_title("Raw signal spectrum")


fig2, axs2 = plt.subplots(3, 1, figsize=(7, 15))

# Print butter_filtered signal
sos = butter(ORDER, Wn=FC, btype='lowpass', analog=False, fs=FR, output='sos')
filtered = sosfilt(sos, signal)
axs2[0].plot(t, filtered, label='Filtered signal', color='C3')
axs2[0].set_title("Filtered signal waveform")
axs2[0].grid(color='grey', linewidth=1, linestyle='--',)

# Print dfft spectrum of butter_filtered signal
yf = np.fft.fft(filtered)
xf = np.fft.fftfreq(len(yf), 1/FR)[:N//2]

axs2[1].plot(xf, np.abs(yf[0:N//2]), color='C5')
axs2[1].axis([0, 2000, 0, 17000])
axs2[1].set_title("Filtered signal spectrum (lowpass butter "+ str(ORDER) +"th order - fc:"+str(FC)+"Hz")

# Print butter filtered signal
for FC in [100, 150, 200]:
  b, a = butter(N=ORDER, Wn=FC, btype='lowpass', analog=False, output='ba', fs=FR)
  w, h = freqz(b=b, a=a, fs=FR, worN=512)
  axs2[2].plot(w, abs(h), label="FC = %d" % FC)

axs2[2].axis([0, 2000, 0, 1])
axs2[2].grid(color='grey', linewidth=1, linestyle='--',)

plt.show()