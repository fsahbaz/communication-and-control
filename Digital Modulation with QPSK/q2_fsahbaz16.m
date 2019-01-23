clc;
clear all;
% Digital signal modulation and demodulation with QPSK

% Part a

% Generating a pseudo-random bipolar binary data stream
seq = ltePRBS(123,10,'signed');
seq = seq.';
% Plotting the initial bipolar binary data stream
figure(1)
stem(seq,'LineWidth',2); grid on;
xlabel('n');
title('Bipolar Bit Sequence');

% Part b

% Setting further parameters
f=1e4; % Bit rate
T=1/f;
% Determining the In-phase and Quadrature components of the modulated
% signal
IC=[];
QS=[];
I=[];
Q=[];
t = [];
tPer = 0:T/100:T-T/100;
for i=1:2:length(seq)  
    if seq(i)==-1 && seq(i+1)==-1
        % Phase pi/4 corresponds to bit sequence -1-1
        I=[I cos(pi/4).*ones(1,length(tPer))]; % I equals 1*cos[phi[n]]
        Q=[Q sin(pi/4).*ones(1,length(tPer))]; % Q equals 1*sin[phi[n]]
        % For the direct calculation of s(t).
        IC=[IC cos(pi/4).*cos(2*pi*f*tPer)]; % I multiplied by cos(2*pi*f*t)
        QS=[QS sin(pi/4).*sin(2*pi*f*tPer)]; % Q multiplied by sin(2*pi*f*t)
    elseif seq(i)==-1 && seq(i+1)==1
        % Phase 3*pi/4 corresponds to bit sequence -11
        Q=[Q sin(3*pi/4).*ones(1,length(tPer))]; % Q equals 1*sin[phi[n]]
        I=[I cos(3*pi/4).*ones(1,length(tPer))]; % I equals 1*cos[phi[n]]
        % For the direct calculation of s(t).
        QS=[QS sin(3*pi/4).*sin(2*pi*f*tPer)]; % Q multiplied by sin(2*pi*f*t)
        IC=[IC cos(3*pi/4).*cos(2*pi*f*tPer)]; % I multiplied by cos(2*pi*f*t)
    elseif seq(i)==1 && seq(i+1)==-1
        % Phase 5*pi/4 corresponds to bit sequence 1-1
        I=[I cos(5*pi/4).*ones(1,length(tPer))]; % I equals 1*cos[phi[n]]
        Q=[Q sin(5*pi/4).*ones(1,length(tPer))]; % Q equals 1*sin[phi[n]]
        % For the direct calculation of s(t).
        IC=[IC cos(5*pi/4).*cos(2*pi*f*tPer)]; % I multiplied by cos(2*pi*f*t)
        QS=[QS sin(5*pi/4).*sin(2*pi*f*tPer)]; % Q multiplied by sin(2*pi*f*t)
    elseif seq(i)==1 && seq(i+1)==1
        % Phase 7*pi/4 corresponds to bit sequence 11
        Q=[Q sin(7*pi/4).*ones(1,length(tPer))]; % Q equals 1*sin[phi[n]]
        I=[I cos(7*pi/4).*ones(1,length(tPer))]; % I equals 1*cos[phi[n]]
        % For the direct calculation of s(t).
        QS=[QS sin(7*pi/4).*sin(2*pi*f*tPer)]; % Q multiplied by sin(2*pi*f*t)
        IC=[IC cos(7*pi/4).*cos(2*pi*f*tPer)]; % I multiplied by cos(2*pi*f*t)
    end
    t = [t tPer];
    tPer = tPer+2;
end
% The modulated output equals I(t)*cos(2*pi*f*t)-Q(t)*sin(2*pi*f*t)
sMod = IC-QS;
time = 0:T:length(sMod)*T-T;
time = time/50;
% Plotting the modulated signal in time domain.
figure(2)
subplot(3,1,1)
plot(time,sMod,'LineWidth',2); grid on;
xlabel('Time (s)');
ylabel('s(t)');
title('QPSK Modulated Message');
% Plotting the in-phase component of the modulated signal in time domain.
subplot(3,1,2)
plot(time,IC); grid on;
hold on;
plot(time,I,'r','LineWidth',2);
xlabel('Time (s)');
ylabel('I(t)');
title('In-phase Component of the Modulated Signal');
% Plotting the quadrature component of the modulated signal in time domain.
subplot(3,1,3)
plot(time,QS); grid on;
hold on;
plot(time,Q,'r','LineWidth',2);
xlabel('Time (s)');
ylabel('Q(t)');
title('Quadrature Component of the Modulated Signal');
% Testing out the output with I and Q components
% figure(3)
% plot(time,I.*cos(2*pi*f*t)-Q.*sin(2*pi*f*t),'LineWidth',2); grid on;
% xlabel('Time (s)');
% ylabel('s(t)');
% title('I(t)cos(\Omega_ct)-Q(t)sin(\Omega_ct)');

% Part c
% Resetting the periodic time interval value in order to use it in periodic integration
tPer = 0:T/100:T-T/100;
mDemod=[];
ICCI=[];
ICC=[];
QSS=[];
for i=1:1:length(seq)/2
    % Multiplying the signal with cosine and integrating to separate the I component
    ICCI = trapz(tPer,sMod((i-1)*length(tPer)+1:i*length(tPer)).*cos(2*pi*f*tPer));
    % Multiplying the signal with sine and integrating to separate the Q component
    QSSI = trapz(tPer,sMod((i-1)*length(tPer)+1:i*length(tPer)).*sin(2*pi*f*tPer));
    % Generating the filter
    sinPart = sin(pi*tPer)./(pi*tPer);
    sinPart(abs(pi*tPer)==0) = 1;
    % Generating the cos and variable parts for Beta=0
    beta = 0;
    cosPart = cos(beta*pi*tPer)./(1-(2*beta*tPer).^2);
    % Multiplying two parts
    g = sinPart.*cosPart;
    gF = fft(g);
    % Passing the signal components through the filter
    ICCIF = fft(ICCI);
    QSSIF = fft(QSSI);
    ICCI = ifft(gF.*ICCIF);
    QSSI = ifft(gF.*QSSIF);
    % Comparing the intgrated and filtered signal components with 0 to further separate it.
    if ICCI>0 
        ICC=[ICC (1/(2^(1/2))*ones(1,length(tPer)))];
    else
        ICC=[ICC (-1*1/(2^(1/2)))*ones(1,length(tPer))];
    end
    
    if QSSI>0
        QSS=[QSS (-1*1/(2^(1/2)))*ones(1,length(tPer))];
    else
        QSS=[QSS (1/(2^(1/2)))*ones(1,length(tPer))];
    end
end
% Obtaining the original bit information by first multiplying the recovered
% I component by two (since both phase pi/d and phase 5*pi/4 components 
% yield 0 when subtracted), and then subtracting Q from it. This way, each
% bit sequence has a corresponding subtracted value to be obtained from it.
for i=1:length(seq)/2
    
    dCON=2.*ICC-QSS;
    if(dCON(i*length(tPer)-length(tPer)/2))==(1/(2^(1/2)))
        mDemod=[mDemod [-1 -1]];
    elseif (dCON(i*length(tPer)-length(tPer)/2))==(-3/(2^(1/2)))
        mDemod=[mDemod [-1 1]];
    elseif (dCON(i*length(tPer)-length(tPer)/2))==(-1/(2^(1/2)))
        mDemod=[mDemod [1 -1]];
    elseif  (dCON(i*length(tPer)-length(tPer)/2))==(3/(2^(1/2)))
        mDemod=[mDemod [1 1]];
    end
    
end
figure(4)
subplot(3,1,1)
plot(time,ICC,'LineWidth',2); grid on;
xlabel('Time (s)');
ylabel('I(t)');
title('Recovered In-Phase Component');
subplot(3,1,2)
plot(time,QSS,'LineWidth',2); grid on;
xlabel('Time (s)');
ylabel('Q(t)');
title('Recovered Q(t) Component');
subplot(3,1,3)
stem(mDemod,'LineWidth',2); grid on;
xlabel('n');
title('Bit Sequence Recovered from the Demodulated Signal');
    