% EECE 4572 Communication Systems
% Homework 2
% Tim Liming

% Variables
f1 = 200; % Hz
f2 = 2000; % Hz
A1 = 1;
A2 = 0.1;
fc = 10000; % Hz
fs = 4*fc; % Hz
dur = 0.01; % s

t = (0:dur*fs)/fs; % time vector

%% 1. Generate u(t) & s(t)

u_t = A1*cos(2*pi*f1*t)+A2*cos(2*pi*f2*t);
s_t = u_t.*cos(2*pi*fc*t);

%% 2. Plot u(t) & s(t)

% Plot of u(t)
figure(1);
plot(t,u_t);
title('u(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot of s(t)
figure(2);
plot(t,s_t);
title('s(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3. Generate U(f) & S(f)

NFFT = 2000; % N-point fft
f = fs/2*linspace(0,1,NFFT/2); % create frequency vector
new_f = cat(2,-1*fliplr(f),f); % invert and flip frequency vector and concatenate
U_f = fftshift(fft(u_t,NFFT)); % U(f)
S_f = fftshift(fft(s_t,NFFT)); % S(f)

%% 4. Plot the spectra

% Plot of U(f)
figure(3);
plot(new_f,abs(U_f));
title('Spectrum of u(t)');
xlabel('Frequency (Hz)');
ylabel('|U(f)|');

% Plot of S(f)
figure(4);
plot(new_f,abs(S_f));
title('Spectrum of s(t)');
xlabel('Frequency (Hz)');
ylabel('|S(f)|');

%% 5. Sounds

ts=(0:fs)/fs; % need to increase the duration
u_t= A1*cos(2*pi*f1*ts)+A2*cos(2*pi*f2*ts); % recalculate u(t)
s_t= u_t.*cos(2*pi*fc*ts); % recalculate s(t)

sound(u_t, fs); % play u(t)
sound(s_t, fs); % play s(t)

%% 6. Experiment

% add a third tone
f3 = 800; % Hz
A3 = 0.8;

ut_ex = A1*cos(2*pi*f1*t)+A2*cos(2*pi*f2*t)+A3*cos(2*pi*f3*t);
st_ex = ut_ex.*cos(2*pi*fc*t);

% Plot of new u(t)
figure(5);
plot(t,ut_ex);
title('Experimental u(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot of new s(t)
figure(6);
plot(t,st_ex);
title('Experimental s(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');

% increase A2
A2 = 0.6;

ut_ex = A1*cos(2*pi*f1*t)+A2*cos(2*pi*f2*t)+A3*cos(2*pi*f3*t);
st_ex = ut_ex.*cos(2*pi*fc*t);

% Plot of new u(t)
figure(7);
plot(t,ut_ex);
title('Increased A2 Experimental u(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot of new s(t)
figure(8);
plot(t,st_ex);
title('Increased A2 Experimental s(t) vs. Time');
xlabel('Time (s)');
ylabel('Amplitude');