% Analog signal amplitude modulation using MATLAB

clc
clear all

% Part (a)

% Reading the audio signal.
[m,fs] = audioread('speech_dft_8kHz.wav');
fc = 160000;
ts = 1/fs;
% Setting the time domain.
time = 0:ts:(length(m)*ts)-ts;
% Upscaling the message signal, and resetting the time domain.
mUp = interp(m,100);
timeUp = 0:ts:(length(mUp)*ts)-ts;
% Modulating the message signal.
y = (2+mUp.').*(cos(2*pi*fc*timeUp));
% Taking the Fourier transform of the signal.
yF = fft(y);
% Setting the frequency domain.
f = fs*(0:length(mUp)-1)/length(mUp);

figure(1)
% Shifting zero-frequency component to center while plotting.
plot(f,mag2db(abs(fftshift(yF))))
xlabel('Frequency (Hz)');
ylabel('|Y(j\Omega)| in dB');
title 'Modulated Message Signal in Fourier Domain';

% Part (b)

tau = 1e-4;
h = exp(-timeUp/tau);
hF = fft(h);
% Passing the differentiated signal through a diode.
y(y<0) = 0;
% Now that the negative parts of the signal is eliminated, it can pass
% through the filter. In order to pass the signal through the filter, both
% of the signal and the transfer function can be multiplied in the
% frequency domain (Fourier transformed).
yZF = fft(y);
yEnv = yZF.*hF;
% Also transforming the message signal in order to compare it with the
% envelope detected output.
mF = fft(mUp);

figure(2)
subplot(2,1,1)
% Shifting zero-frequency component to center while plotting.
plot(f,mag2db(abs(fftshift(yEnv))))
xlabel('Frequency (Hz)');
ylabel('|Mj\Omega)| (Demodulated) in dB');
title 'Demodulated Message Signal in Fourier Domain';
subplot(2,1,2)
% Shifting zero-frequency component to center while plotting.
plot(f,mag2db(abs(fftshift(mF))))
xlabel('Frequency (Hz)');
ylabel('|M(j\Omega)| (Original) in dB');
title 'Message Signal in Fourier Domain';

mD = ifft(yEnv);

figure(3)
subplot(2,1,1)
plot(timeUp,mD);
xlabel('Time(s)');
ylabel('|m(t)| (Demodulated)');
title 'Demodulated Message Signal (Not downsampled yet)';
subplot(2,1,2)
plot(time,m);
xlabel('Time (s)');
ylabel('|m(t)| (Original)');
title 'Message Signal';

% Downsampling the signal by 100.
mDD = decimate(mD,100);
% Subtracting 2 in order to recover the original signal.

mDD = mDD-2;
sound(mDD);

% Part (c)

% Modulating the message signal.
yS = (mUp.').*(cos(2*pi*fc*timeUp));
% Taking the Fourier transform of the modulated signal.
ySF = fft(yS);
% Passing the differentiated signal through a diode.
yS(yS<0) = 0;
% Now that the negative parts of the signal is eliminated, it can pass
% through the filter. In order to pass the signal through the filter, both
% of the signal and the transfer function can be multiplied in the
% frequency domain (Fourier transformed).
ySZF = fft(yS);
ySEnv = ySZF.*hF;
% Taking the inverse Fourier transform of the envelope detected output.
mSD = ifft(ySEnv);

figure(4)
subplot(2,1,1)
% Shifting zero-frequency component to center while plotting.
plot(timeUp,mSD)
xlabel('Time (s)');
ylabel('|m(t)| (Demodulated)');
title 'Demodulated Message Signal (Not downsampled yet)';
subplot(2,1,2)
% Shifting zero-frequency component to center while plotting.
plot(time,m)
xlabel('Frequency (Hz)');
ylabel('|m(t)| (Original)');
title 'Message Signal';

% Downsampling the signal by 100.
mSDD = decimate(mSD,100);
sound(mSDD);

