% EECE4572 Communication Systems
% Tim Liming
% Homework 7

%% Variables

num_bits = 1000;    % number of bits in the bit sequence
bit_seq = round(rand([1 num_bits]));	% logical bit sequence
R = 1000;       % symbols per second
T = 1/R;        % pulse duration
N = 16;         % number of samples
fs = N*R;       % sampling frequency
fc = 100*R;     % carrier frequency
Eg = 1;         % energy
A = sqrt(Eg/T); % amplitude
t = (0:1/fs:(R/2)*T);   % time vector

theta_e = 0;    % theta error

%% 4-QAM

% build d(n)
dn1 = zeros([1 round(numel(bit_seq)/2)]);
j = 1;  % initialize index var
for i = 1:2:numel(bit_seq)
    if     bit_seq(i:i+1) == [0,0]
        dn1(j) = 1+1i;
    elseif bit_seq(i:i+1) == [0,1]
        dn1(j) = -1+1i;
    elseif bit_seq(i:i+1) == [1,1]
        dn1(j) = -1-1i;
    elseif bit_seq(i:i+1) == [1,0]
        dn1(j) = 1-1i;
    end % if
    j = j+1;
end % for

% generate baseband signal u(t)
ut1 = [];           % initialize u(t)
gt = A*ones(1,N);   % create g(t)
for i = 1:numel(dn1)
    ut1 = cat(2,ut1,dn1(i)*gt);
end % for

% plot u(t)
figure(1);
%plot(t(1:end-1),ut1);
plot(ut1,'*');
title('Baseband Signal u(t)  (4-QAM)');

% generate modulated signal s(t)
st1 = real(ut1.*exp(1i*2*pi*fc*t(1:end-1)));

% calculate r(t)
temp = (st1.*cos(2*pi*fc*t(1:end-1)))*2;
rt1 = [];    % initialize r(t)
for i = 1:N:numel(temp)
    rt1 = [rt1 mean(temp(i:i+(N-1)))];
end % for

% plot recieved signal
figure(2);
plot(rt1,'*');
title('Recovered signal r(t)  (4-QAM)');
xlabel('Symbols');
ylabel('r(t)');

% check recieved signal
figure(3);
plot(dn1,rt1);
title('Recovered Symbols vs. Symbols  (4-QAM)');

%% QPSK (4-PSK)

% build d(n)
dn2 = zeros([1 round(numel(bit_seq)/2)]);
j = 1;  % initialize index var
for i = 1:2:numel(bit_seq)
    if     bit_seq(i:i+1) == [0,0]
        dn2(j) = exp(1i*0);
    elseif bit_seq(i:i+1) == [0,1]
        dn2(j) = exp(1i*pi/2);
    elseif bit_seq(i:i+1) == [1,1]
        dn2(j) = exp(1i*pi);
    elseif bit_seq(i:i+1) == [1,0]
        dn2(j) = exp(1i*3*pi/2);
    end % if
    j = j+1;
end % for

% generate baseband signal u(t)
ut2 = [];           % initialize u(t)
gt = A*ones(1,N);   % create g(t)
for i = 1:numel(dn2)
    ut2 = cat(2,ut2,dn2(i)*gt);
end % for

% plot u(t)
figure(4);
%plot(t(1:end-1),ut2);
plot(ut2,'*');
title('Baseband Signal u(t)  (QPSK)');

% generate modulated signal s(t)
st2 = real(ut2.*exp(1i*2*pi*fc*t(1:end-1)));

% calculate r(t)
temp = (st2.*cos(2*pi*fc*t(1:end-1)))*2;
rt2 = [];    % initialize r(t)
for i = 1:N:numel(temp)
    rt2 = [rt2 mean(temp(i:i+(N-1)))];
end % for

% plot recieved signal
figure(5);
plot(rt2,'*');
title('Recovered signal r(t)  (QPSK)');
xlabel('Symbols');
ylabel('r(t)');

% check recieved signal
figure(6);
plot(dn2,rt2);
title('Recovered Symbols vs. Symbols  (QPSK)');

%% Differential QPSK  (DQPSK)

% build d(n)
dn3 = zeros([1 round(numel(bit_seq)/2)]);
j = 1;  % initialize index var
for i = 1:2:numel(bit_seq)
    if     bit_seq(i:i+1) == [0,0]
        dn3(j) = exp(1i*0);
    elseif bit_seq(i:i+1) == [0,1]
        dn3(j) = exp(1i*pi/2);
    elseif bit_seq(i:i+1) == [1,1]
        dn3(j) = exp(1i*pi);
    elseif bit_seq(i:i+1) == [1,0]
        dn3(j) = exp(1i*3*pi/2);
    end % if
    j = j+1;
end % for

% generate baseband signal u(t)
ut3 = [];           % initialize u(t)
gt = A*ones(1,N);   % create g(t)
for i = 1:numel(dn3)
    ut3 = cat(2,ut3,dn3(i)*gt);
end % for

% plot u(t)
figure(7);
%plot(t(1:end-1),ut2);
plot(ut3,'*');
title('Baseband Signal u(t)  (DQPSK)');

% generate modulated signal s(t)
st3 = real(ut3.*exp(1i*2*pi*fc*t(1:end-1)));

% calculate r(t)
temp = (st3.*cos(2*pi*fc*t(1:end-1)))*2;
rt3 = [];    % initialize r(t)
for i = 1:N:numel(temp)
    rt3 = [rt3 mean(temp(i:i+(N-1)))];
end % for

% plot recieved signal
figure(8);
plot(rt1,'*');
title('Recovered signal r(t)  (DQPSK)');
xlabel('Symbols');
ylabel('r(t)');

% check recieved signal
figure(9);
plot(dn3,rt3);
title('Recovered Symbols vs. Symbols  (DQPSK)');
