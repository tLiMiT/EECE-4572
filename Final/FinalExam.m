% EECE4572 Communication Systems
% Final Exam
% Tim Liming

%% Problem 1

% BPSK system
load('problem1_13.mat');
T = 1*10^-3; % T = 1ms
R = 16; % sample rate
fs = R*10^3; % 16kHz
t = 10^3*((1/fs):(1/fs):length(v)/fs); % time vector

% input signal plot
figure(1);
plot(t,v);
title('problem1\_13.mat');
xlabel('time (ms)');
ylabel('noisy signal');

% matched filter
%mf = conv(v, fliplr(conj((v))), 'same');

% recover the bits
rb = [];
index = 1;
for i = 1:R:length(v)
    % This is our "matched filter"
    rb(index) = mean(v(i:i+(R-1)));
    index = index +1;
end % for

% build bit string
bs = [];
for i = 1:length(rb)
    if (rb(i) < 0)
        bs(i) = 0;
    else
        bs(i) = 1;
    end % if
end % for

% match symbols
msg = [];
for i = 1:5:length(bs)
    % assign letters
    switch (bin2dec(num2str(bs(i:i+4))))
        % bits have been translated to decimal values
        case 0
            msg = strcat(msg,{' '});
        case 1
            msg = strcat(msg,'a');
        case 2
            msg = strcat(msg,'b');
        case 3
            msg = strcat(msg,'c');
        case 4
            msg = strcat(msg,'d');
        case 5
            msg = strcat(msg,'e');
        case 6
            msg = strcat(msg,'f');
        case 7
            msg = strcat(msg,'g');
        case 8
            msg = strcat(msg,'h');
        case 9
            msg = strcat(msg,'i');
        case 10
            msg = strcat(msg,'j');
        case 11
            msg = strcat(msg,'k');
        case 12
            msg = strcat(msg,'l');
        case 13
            msg = strcat(msg,'m');
        case 14
            msg = strcat(msg,'n');
        case 15
            msg = strcat(msg,'o');
        case 16
            msg = strcat(msg,'p');
        case 17
            msg = strcat(msg,'q');
        case 18
            msg = strcat(msg,'r');
        case 19
            msg = strcat(msg,'s');
        case 20
            msg = strcat(msg,'t');
        case 21
            msg = strcat(msg,'u');
        case 22
            msg = strcat(msg,'v');
        case 23
            msg = strcat(msg,'w');
        case 24
            msg = strcat(msg,'x');
        case 25
            msg = strcat(msg,'y');
        case 26
            msg = strcat(msg,'z');
    end % switch-case
end % for

display(msg);

%% Problem 2

% BPSK system
Eb_No = 8; % SNR (dB)
theta_e = pi; % phase offset
ang_range = -theta_e:(pi/16):theta_e;
% probability of error
Pe = qfunc(sqrt(2*db2pow(Eb_No))*real(exp(1i*ang_range)));
% DBPSK system
Pb = (1/2)*exp(-1*db2pow(Eb_No)); % bit error
invPb = qfuncinv(Pb);
angle = invPb/sqrt(2*db2pow(Eb_No));

P2ans = (acos(angle)); % radians
display(P2ans);

%% Problem 3

%QAM system
M = 4^4; % M-ary
SNR = 23; % dB (Eb/No)
f = 1*10^6; % 1 MHz
Tinv = f/2; % B=2/T -> 1/T
BER = 10^-5; % BER (without coding)

% bit error
Pb = (1/log2(M))*4*((sqrt(M)-1)/(sqrt(M)))*...
    qfunc(sqrt(3*log2(M)/(M-1)*db2pow(SNR)));

% bit rate
Rb = (log2(M)*Tinv)/10^6; % Mbits
display(Rb); % Mbits

% less than BER?
less = Pb < BER; % true (1) or false (0)
display(less);

%% Problem 4

% BFSK
Rb = 25*10^6; % bps
Pt = 200; % W (transmission power)
A = 140; % dB (attenuation)
NoE = 4*10^-21; % W/Hz (@earth)

% a) regenerator
% noise at satellite reciever is 10dB lower
No = NoE/db2pow(10); % No (@space)

% at the satellite (point a)
Pna = Rb*No;
PRa = Pt*db2pow(-A)/Pna;
Pea = (1/2)*exp(-1*PRa/2);

% back on earth (point b)
Pnb = Rb*NoE;
PRb = Pt*db2pow(-A)/Pnb;
Peb = (1/2)*exp(-1*PRb/2);

% Total probability error for regenerator
PeA = Pea + Peb;
display(PeA);


% b) repeater
% back on earth (point b)
PRb = Pt*db2pow(-A)/(Pnb+Pna);
Peb = (1/2)*exp(-1*PRb/2);

% Total probability error for repeater
PeB = Pea + Peb;
display(PeB);
