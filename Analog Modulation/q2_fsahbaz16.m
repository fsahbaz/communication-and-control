% Analog signal frequency modulation using MATLAB

clc
clear all

% Part (a)

% Reading the audio signal.
[m,fs] = audioread('speech_dft_8kHz.wav');
fc = 280000;
kf = 75000;
ts = 1/fs;
% Setting the time domain.
time = 0:ts:(length(m)*ts)-ts;
% Upsampling the message signal, and resetting the time domain.
mUp = interp(m,200);
timeUp = 0:ts:(length(mUp)*ts)-ts;
% Modulating the message signal.
y = cos(2*pi*fc*timeUp + (ts*kf/200)*cumsum(mUp.')) ;
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

tau = 1e-3;
h = exp(-timeUp/tau);
hF = fft(h);
% Differentiating the frequency modulated signal and multiplying it by the 
% sampling rate in order to rescale it.
yD = diff(y)*(200*fs/kf);
% Passing the differentiated signal through a diode.
yD(yD<0) = 0;
% Now that the negative parts of the signal is eliminated, it can pass
% through the filter. In order to pass the signal through the filter, both
% of the signal and the transfer function can be multiplied in the
% frequency domain (Fourier transformed).
yDZF = fft(yD);
yEnv = hF(1:end-1).*yDZF;
% Also transforming the message signal in order to compare it with the
% envelope detected output.
mF = fft(mUp);

figure(2)
subplot(2,1,1)
% Shifting zero-frequency component to center while plotting.
plot(f(1:end-1),mag2db(abs(fftshift(yEnv))))
xlabel('Frequency (Hz)');
ylabel('|M(j\Omega)| (Demodulated) in dB');
title 'Demodulated Message Signal in Fourier Domain';
subplot(2,1,2)
% Shifting zero-frequency component to center while plotting.
plot(f,mag2db(abs(fftshift(mF))))
xlabel('Frequency (Hz)');
ylabel('|M(j\Omega)| (Original) in dB');
title 'Message Signal in Fourier Domain';

% Taking the inverse Fourier transform of the envelope detected output.
mD = ifft(yEnv)./8;


figure(3)
subplot(2,1,1)
plot(timeUp(1:end-1),abs(mD));
xlabel('Time(s)');
ylabel('|m(t)| (Demodulated)');
title 'Demodulated Message Signal';
subplot(2,1,2)
plot(time,m);
xlabel('Time (s)');
ylabel('|m(t)| (Original)');
title 'Message Signal';

% Downsampling the signal by 200.
mDD = decimate(mD,200);
soundsc(real(mDD));