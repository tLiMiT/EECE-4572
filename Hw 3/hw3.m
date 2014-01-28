% EECE4572 Communication Systems
% Tim Liming
% Homework 3

% Variables
A = 1; % amplitude
fc = 10E3; % Hz (carrier)
fm = fc/1000; % Hz (modulator)
fs = 4*fc; % Hz (sampling)
delta_t = 1/fs; % time step
T = 1000/fc; % time interval
t = (0:delta_t:T); % time
Pm = 3; % phi-m 

%% Input

P = Pm*cos(2*pi*fm*t); % phi(t) (eqn 1)
r = A*cos(2*pi*fc*t+P); % r(t) (input to PLL)

%% PLL

% Variables
M = 2*fs/fc; 
K1 = 0.1; % PLL Constant 1
K2 = K1/10; % PLL Constant 2

% Set initial values to 0
P_hat = zeros(1,length(t)); 
Pe_tilde = zeros(1,length(t)); 
Pe = zeros(1,length(t));

% loop for PLL
for n=1:length(t)
    Pe_tilde(n) = -2*sin(2*pi*fc*t(n)+P_hat(n)).*r(n); % (eqn 2)
    Pe(n) = (1/M)*sum(Pe_tilde(max(n-M+1,1):n)); % Low Pass Filter (eqn 3)
    P_hat(n+1) = P_hat(n)+K1*Pe(n)+K2*sum(Pe); % (eqn 4)
end % for

%% Plots

figure(1);
plot(t,P,'-r', t,P_hat(1:length(P_hat)-1),'-b');
title('\phi(n\Deltat), \phi(n\Deltat)-hat vs Time');
xlabel('Time (s)');
ylabel('Amplitude');
legend('\phi(n\Deltat)', '\phi(n\Deltat)-hat');