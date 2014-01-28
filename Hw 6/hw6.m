% EECE4572 Communication Systems
% Tim Liming
% Homework 6

%% Variables

bit_seq = round(rand([1 1000]));    % logical bit sequence
R = 1000;           % symbol rate (symbols per second)
T = 1/R;            % pulse duration
N = 400;            % number of samples
fs = N*R;           % sampling frequency (Hz)
Ts = 1/fs;          % Sample period
fc = 100*R;         % carrier frequency
Eg = 1;             % energy
A = sqrt(Eg/T);     % amplitude
t = (0:1/fs:R*T);   % time vector

%% Binary PAM

M1 = 2;

% build d(n)
dn2 = zeros([1 numel(bit_seq)]);
for i = 1:numel(bit_seq)
    if bit_seq(i) == 1
        dn2(i) = 1;
    else
        dn2(i) = -1;
    end % if
end % for

% calculate u(t)
ut = [];            % initialize u(t)
gt = A*ones(1,N);   % create g(t)
for i = 1:numel(dn2)
    ut = cat(2,ut,dn2(i)*gt);
end % for

% plot u(t)
figure(1);
plot(t(1:end-1),ut);
title('Baseband Signal u(t)  (M = 2)');
xlabel('Time (s)');
ylabel('Amplitude');

% calculate the signal spectrum
f = fs/2*linspace(0,1,fs/2);    % create frequency vector
new_f = cat(2,-1*fliplr(f),f);  % invert and flip frequency vector and cat
U_f = fftshift(fft(ut,fs));     % U(f)

% Plot of Su(f)
figure(2);
plot(new_f,abs(U_f));
title('Spectrum of u(t)  (M = 2)');
xlabel('Frequency (Hz)');
ylabel('Su(f)');

% calculate s(t)
st = ut.*cos(2*pi*fc*t(1:end-1));

% calculate the signal spectrum
S_f = fftshift(fft(st,fs));  % S(f)

% Plot of Ss(f)
figure(3);
plot(new_f,abs(S_f));
title('Spectrum of s(t)  (M = 2)');
xlabel('Frequency (Hz)');
ylabel('Ss(f)');

% calculate r(t)
temp = (st.*cos(2*pi*fc*t(1:end-1)))*2;
rt = [];    % initialize r(t)
for i = 1:N:numel(temp)
    rt = [rt mean(temp(i:i+(N-1)))];
end % for

% plot recieved signal
figure(4);
plot(rt);
title('Recovered signal r(t)  (M = 2)');
xlabel('Symbols');
ylabel('r(t)');

% check recieved signal
figure(5);
plot(dn2,rt);
title('Recovered Symbols vs. Symbols');


%% M-PAM

M2 = 8;

% build d(n)  (using Gray mapping)
dn8 = zeros([1 round(numel(bit_seq)/3)]);
for i = 1:3:(numel(bit_seq)-1)
    if     bit_seq(i:i+2) == [0,0,0]
        dn8(ceil(i/3)) = 7;         % set to +7
    elseif bit_seq(i:i+2) == [0,0,1]
        dn8(ceil(i/3)) = 5;         % set to +5
    elseif bit_seq(i:i+2) == [0,1,1]
        dn8(ceil(i/3)) = 3;         % set to +3
    elseif bit_seq(i:i+2) == [0,1,0]
        dn8(ceil(i/3)) = 1;         % set to +1
    elseif bit_seq(i:i+2) == [1,1,0]
        dn8(ceil(i/3)) = -1;        % set to -1
    elseif bit_seq(i:i+2) == [1,1,1]
        dn8(ceil(i/3)) = -3;        % set to -3
    elseif bit_seq(i:i+2) == [1,0,1]
        dn8(ceil(i/3)) = -5;        % set to -5
    elseif bit_seq(i:i+2) == [1,0,0]
        dn8(ceil(i/3)) = -7;        % set to -7
    end % if-elseif
end % for

% calculate u(t)
ut2 = [];           % initialize u(t)
gt = A*ones(1,N);   % create g(t)
for i = 1:numel(dn8)
    ut2 = cat(2,ut2,dn8(i)*gt);
end % for

% plot u(t)
t2 = (0:333*N)*Ts;  % time vector
figure(6);
plot(t2(1:end-1),ut2);
title('Baseband Signal u(t)  (M = 8)');
xlabel('Time (s)');
ylabel('Amplitude');

% calculate the signal spectrum
NFFT = numel(t2)-1;
f = fs/2*linspace(0,1,NFFT/2);  % create frequency vector
new_f = cat(2,-1*fliplr(f),f);  % invert and flip frequency vector and cat
U_f2 = fftshift(fft(ut2,NFFT));	% U(f)

% Plot of Su(f)
figure(7);
plot(new_f,abs(U_f2));
title('Spectrum of u(t)  (M = 8)');
xlabel('Frequency (Hz)');
ylabel('Su(f)');

% calculate s(t)
st2 = ut2.*cos(2*pi*fc*t2(1:end-1));

% calculate the signal spectrum
S_f2 = fftshift(fft(st2,numel(t2)-1));  % S(f)

% Plot of Ss(f)
figure(8);
plot(new_f,abs(S_f2));
title('Spectrum of s(t)  (M = 8)');
xlabel('Frequency (Hz)');
ylabel('Ss(f)');

% calculate r(t)
temp = (st2.*cos(2*pi*fc*t2(1:end-1)))*2;
rt2 = [];    % initialize r(t)
for i = 1:N:numel(temp)
    rt2 = [rt2 mean(temp(i:i+(N-1)))];
end % for

% plot recieved signal
figure(9);
plot(rt2);
title('Recovered signal r(t)  (M = 8)');
xlabel('Symbols');
ylabel('r(t)');

% check recieved signal
figure(10);
plot(dn8,rt2);
title('Recovered Symbols vs. Symbols');
