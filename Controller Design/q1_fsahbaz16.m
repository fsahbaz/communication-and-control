clear all
clc

% All figures have been exported with log-scale on the x-axis, in order
% for them to be more observable.

% Definitions:
syms s
Kp = 7.583; 
Ki = 0.2906;     
% Generating the plant tranfer function.
TF = tf(1,[1,3,1]);
% Setting the proportional component of the controller.
P = Kp;
% Setting the integral component of the controller.
I = tf(Ki,[1,0]); % Ki/s
% Generating the open-loop plant.
S = TF; % feedback(TF,1);--will not use this one since the inital plant is
% an open-loop.
% Plotting the step response of the plant.
[yS,t1] = step(S,10000);
figure(1)
% step(S)
subplot(2,1,1)
plot(t1,yS,'LineWidth',2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of the Open-Loop Plant');
% Plotting the root locus representation.
subplot(2,1,2)
rlocus(S)
% Plotting P controller closed-loop step response for various Kp values.
figure(2)
% Kp = 0.4739375
P = 0.4739375;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'--r');
hold on;
% Kp = 0.947875
P = 0.947875;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'--r');
hold on;
% Kp = 1.89575
P = 1.89575;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'--r');
hold on;
% Kp = 3.7915
P = 3.7915;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'--r');
hold on;
% Kp = 7.583
P = 7.583;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'g','LineWidth',2);
hold on;
% Kp = 15.166
P = 15.166;
S_P = feedback(TF*P,1);
[ySP,t2] = step(S_P,10000);
plot(t2,ySP,'--r');
% Titles
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of the Closed-Loop Plant with Proportional Control');
figure(3)
% Final P Closed-Loop w/ P Controller plot.
% Generating the closed-loop plant with P control.
P = Kp;
S_P = feedback(TF*P,1);
% Plotting the step response of the plant.
[ySP,t2] = step(S_P,10000);
% step(S_P);
subplot(2,1,1)
plot(t2,ySP,'LineWidth',2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of the Closed-Loop Plant with Proportional Control');
% Plotting the root locus representation.
subplot(2,1,2)
rlocus(S_P)
% Plotting P controller closed-loop step response for various Ki values.
figure(4)
% Ki = 0.0181625
I = tf(0.0181625,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'--r');
hold on;
% Ki = 0.036325
I = tf(0.036325,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'--r');
% Ki = 0.07265
I = tf(0.07265,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'--r');
% Ki = 0.1453
I = tf(0.1453,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'--r');
% Ki = 0.2906
I = tf(0.2906,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'g','LineWidth',2);
% Ki = 0.5812
I = tf(0.5812,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
[ySPI,t3] = step(S_PI,10000);
plot(t3,ySPI,'--r');
% Titles
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of the Closed-Loop Plant with Proportional-Integral Control');
% Final PI Controller w\ Closed-Loop plot.
% Generating the closed-loop plant with PI control.
I = tf(Ki,[1,0]);
S_PI = feedback(TF*(P+I),1); % 1 for unity feedback
% Plotting the step response of the plant.
[ySPI,t3] = step(S_PI,10000);
figure(5)
% step(S_PI);
subplot(2,1,1)
plot(t3,ySPI,'LineWidth',2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of the Closed-Loop Plant with Proportional-Integral Control');
% Plotting the root locus representation.
subplot(2,1,2)
rlocus(S_PI)
% Extracting information of step responses of the plants.
S1 = stepinfo(S);
S2 = stepinfo(S_P);
S3 = stepinfo(S_PI);
% Calculating steady-state errors.
SE1 = abs(1-yS(end));
SE2 = abs(1-ySP(end));
SE3 = abs(1-ySPI(end));
% Converting the continuous time controller to a digital controller.
SD = c2d((P+I),0.001,'tustin'); % tustin (bilinear approximation) method.
% Rebuilding the system with a D/A converter
% Also introducing a delay to the plant, in order to realize the digital
% controller and avoid variable delays, since the A/D and D/A converters
% cause delays in real life systems. Padé helps introduce this delay
% without further errors. (Some errors occured when I directly multiplied
% with exp(-0.001*s))
s = tf('s');
TF = tf(1,[1,3,1],'InputDelay',0.001);
PID = feedback(pade(TF)*d2c(SD,'zoh'),1);
[ySPID,t4] = step(PID,1000);
% Plotting the response and the root locus representation of the digital system.
figure(6)
subplot(2,1,1)
plot(t4,ySPID,'--','LineWidth',2);
xlabel('n');
ylabel('Amplitude');
title('Step Response of the Discrete Time Closed-Loop Plant with Proportional-Integral Control');
subplot(2,1,2)
rlocus(PID);
% Extracting information of step responses of the plants.
SD3 = stepinfo(ySPID); % Time related values will be divided by 1000
% since the sampling period was 0.001.
% Calculating steady-state errors.
SDE3 = abs(1-ySPID(end));