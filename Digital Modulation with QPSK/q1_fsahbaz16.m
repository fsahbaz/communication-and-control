% Generating a Raised Cosine Pulse

% Generating the sinc function
fs = 5;
t = -fs*5:1/fs:fs*5;
sinPart = sin(pi*t)./(pi*t);
sinPart(abs(pi*t)==0) = 1;
% Generating the cos and variable parts for Beta=0
beta1 = 0;
cosPart1 = cos(beta1*pi*t)./(1-(2*beta1*t).^2);
% Multiplying two parts
g1 = sinPart.*cosPart1;
gF1 = fft(g1);
f = fs*(0:length(g1)-1)/length(g1);
% Generating the cos and variable parts for Beta=0.5;
beta2 = 0.5;
cosPart2 = cos(beta2*pi*t)./(1-(2*beta2*t).^2);
cosPart2(abs(1-(2*beta2*t).^2)==0)=pi/4;
% Multiplying two parts
g2 = sinPart.*cosPart2;
gF2 = fft(g2);
% Generating the cos and variable parts for Beta=1;
beta3 = 1;
cosPart3 = cos(beta3*pi*t)./(1-(2*beta3*t).^2);
cosPart3(abs(1-(2*beta3*t).^2)==0)=pi/4;
% Multiplying two parts
g3 = sinPart.*cosPart3;
gF3 = fft(g3);

% Part a
% Plotting g in time domain
% Beta=0
figure(1)
plot(t,g1,'LineWidth',1.2)
hold on;
% Beta=0.5
plot(t,g2,'LineWidth',1.2)
hold on;
% Beta=1
plot(t,g3,'LineWidth',1.2)
xlabel('Time (s)');
ylabel('g(t)');
grid on;

% Part b
% Plotting G in frequency domain
% Beta = 0
figure(2)
plot(f,fftshift(abs(gF1)),'LineWidth',1.2)
hold on;
% Beta = 0.5
plot(f,fftshift(abs(gF2)),'LineWidth',1.2)
hold on;
% Beta = 1
plot(f,fftshift(abs(gF3)),'LineWidth',1.2)
xlabel('Frequency (Hz)');
ylabel('G(j\Omega)');
grid on;