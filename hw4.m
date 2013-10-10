% EECE4572 Communication Systems
% Tim Liming
% Homework 4

%% Variables
K = 111;                    % Number of components
fk = (100:30:3400);         % Frequency components
Ak = rand([1 K]);           % Amplitudes
phi_k = rand([1 K])*2*pi;   % Phases
T = 2;                      % signal duration (2s)
fs = 8E3;                   % 8 kHz sampling freq.
Fs = 10*fs;                 % sampling freq. mimicking CT

t = (0:1/Fs:T);             % time vector

%% 1. Generate x(t)

x_t = 0; % Initialize x(t)

for i=1:K
    % Build Signal Vector (Equation 1)
    x_t = x_t + Ak(i)*cos(2*pi*fk(i).*t+phi_k(i));
end % for

% Plot x(t)
figure(1);
plot(t,x_t);
title('x(t) vs Time');
xlabel('Time (s)');
ylabel('x(t)');

%% 2. Sample the signal

I = Fs/fs; % Sample reduction interval

xs = x_t(1:I:end);  % Sampled signal vector
ts = t(1:I:end);    % Sampled time vector

figure(2);
plot(ts, xs);
title('Sampled X vs Time');
xlabel('Time (s)');
ylabel('Sampled X');

%% 3. Quantize the samples

R_N = [3 4 5 6 7 8];        % Quantization levels
N = 2.^R_N;                 % Uniform quantizer levels

xs_normal = xs/(max(abs(xs))*1.000001); % Normalized Xs

for i=1:length(N)
    
    delta = 2/N(i); % delta
    
    % The following generates a vector with the 
    % Quantization levels
%     A = linspace(-1,0,(N(i)/2)+1);
%     B = linspace(0,1,(N(i)/2)+1);
%     Q_levels = cat(2,A(1:end-1),B);
    
    % Quantized samples
    temp = round(xs_normal/delta-0.5);
    q(i,:) = ((temp)+0.5)*delta;

    % Plot the quantized samples
    figure(2+i);
    plot(ts,q(i,:));
    title('Quantized Samples vs Time');
    xlabel('Time (s)');
    ylabel('q');
    
    % Calculate SQNR
    xs_power = mean(power(xs_normal,2));
    xs_q_power = mean(power(xs_normal - q(i,:),2));
    SQNR = xs_power/xs_q_power; % calculate SQNR
    SQNR_dB(i) = pow2db(SQNR);  % convert into dB
    
end % for

%% 4. Plot SQNR

% Plot SQNR
figure(9);
plot(R_N,SQNR_dB, R_N,6*R_N);
title('SQNR and 6*RN vs RN');
xlabel('RN');
ylabel('SQNR');
legend('SQNR','6*RN','Location','NorthWest');

%% 6. Calculate Rb

Rb = R_N*fs;
