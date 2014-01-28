% EECE4572 Communication Systems
% Tim Liming
% Homework 8

%% Problem 2

EbNo = (0:0.1:100);         % bit SNR
Pe = qfunc(sqrt(2*EbNo));	% BER

% plot
figure(1);
semilogy(pow2db(EbNo),(Pe));
title('Pe vs. Eb/No');
xlabel('SNR');
ylabel('BER');

%% i
EbNo1 = 4; % dB
Pe1 = qfunc(sqrt(2*db2pow(EbNo1)));	% BER1
display(Pe1);

% plot
figure(2);
semilogy(pow2db(EbNo),(Pe));
hold on;
plot(EbNo1,Pe1,'-*');
title('Pe vs. Eb/No');
xlabel('SNR');
ylabel('BER')

%% ii
EbNo2 = 8; % dB
Pe2 = qfunc(sqrt(2*db2pow(EbNo2)));	% BER2
display(Pe2);

% plot
figure(3);
semilogy(pow2db(EbNo),(Pe));
hold on;
plot(EbNo2,Pe2,'-*');
title('Pe vs. Eb/No');
xlabel('SNR');
ylabel('BER')

%% iii
EbNo3 = 10; % dB
Pe3 = qfunc(sqrt(2*db2pow(EbNo3)));	% BER3
display(Pe3);

% plot
figure(4);
semilogy(pow2db(EbNo),(Pe));
hold on;
plot(EbNo3,Pe3,'-*');
title('Pe vs. Eb/No');
xlabel('SNR');
ylabel('BER')

%% iv
BER = 10^-5;
SNR = (qfuncinv(BER)^2)/2;
display(SNR);

%% v
No = 4*10^-21;      % W/Hz
lnk_atten = 144;    % dB
Eb = db2pow(pow2db(No)+pow2db(SNR)+lnk_atten);
display(Eb);

%% vi
Rb = 1*10^6; % bps
Pt = Eb*Rb;  % transmit power  (W)
display(Pt);

%% vii
Pt = 1; % W (transmit power)
BER = 10^-5;
RbMax = Pt/Eb;   % maximum bit rate
display(RbMax);

