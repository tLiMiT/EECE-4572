% EECE4572 Communication Systems
% Tim Liming
% Homework 9

%% Variables
N = 100000;
bit_seq = round(rand([1 N]));
Eg = 1;
x = randn([1 N]);
y = randn([1 N]);

% build d(n) for M=2
dn2 = zeros([1 round(numel(bit_seq))]);
for i = 1:numel(bit_seq)
    if     bit_seq(i) == [0]
        dn2(i) = exp(1i*0);
    elseif bit_seq(i) == [1]
        dn2(i) = exp(1i*pi);
    end % if
end % for

% build d(n) for M=4
dn4 = zeros([1 round(numel(bit_seq)/2)]);
j = 1;  % initialize index var
for i = 1:2:numel(bit_seq)
    if     bit_seq(i:i+1) == [0,0]
        dn4(j) = exp(1i*0);
    elseif bit_seq(i:i+1) == [0,1]
        dn4(j) = exp(1i*pi/2);
    elseif bit_seq(i:i+1) == [1,1]
        dn4(j) = exp(1i*pi);
    elseif bit_seq(i:i+1) == [1,0]
        dn4(j) = exp(1i*3*pi/2);
    end % if
    j = j+1;
end % for

%% BPSK

M = 2; % binary
Eb = (1/2)*(1/log2(M))*Eg; % bit energy

Pb1 = [];
EbNo1 = [];
BER1 = [];
for SNR = (0:1:10)
    Eb_No = db2pow(SNR); % Eb/No
    No = Eb/Eb_No; % variance (No*Eg)
    sigma = sqrt(No);
    
    Zc = sigma * x; % real part of noise
    Zs = sigma * y; % imaginary part of noise
    Zn = Zc + 1i*Zs; % total noise
    
    % recovered signal
    yn1 = Eg*(dn2) + (Zn);
    
    decodedBit = zeros([1 length(yn1)]);
    for j = 1:length(yn1)
        if yn1(j) > 0
            decodedBit(j) = 0;
        elseif yn1(j) < 0
            decodedBit(j) = 1;
        end % if
    end % for
    
    correct = 0;
    errors = 0;
    for j = 1:length(yn1)
        if decodedBit(j) == bit_seq(j)
            correct = correct + 1;
        elseif decodedBit(j) ~= bit_seq(j)
            errors = errors + 1;
        end % if
    end % for
    
    total = errors + correct;
    BER1 = [BER1 errors/total];
    
    Pb1 = [Pb1 qfunc(sqrt(2*Eb_No))];
    EbNo1 = [EbNo1 pow2db(Eb_No)];
    
end % for

% plots
figure(1);
plot(yn1, '*');
title('Scatterplot of Symbols at SNR=10dB  (BPSK)');
xlabel('Real');
ylabel('Imaginary');

figure(2);
semilogy(EbNo1, BER1, EbNo1, Pb1);
title('BER vs. Eb/No  (BPSK)');
xlabel('Eb/No (dB)');
ylabel('BER (Log Scale)');
legend('Actual', 'Theoretical');

%% QPSK

M = 4; % quad
Eb = (1/2)*(1/log2(M))*Eg; % bit energy

Pb2 = [];
EbNo2 = [];
BER2 = [];
for SNR = (0:1:10)
    Eb_No = db2pow(SNR); % Eb/No
    No = Eb/Eb_No; % variance (No*Eg)
    sigma = sqrt(No);
    
    Zc = sigma * x; % real part of noise
    Zs = sigma * y; % imaginary part of noise
    Zn = Zc + 1i*Zs; % total noise
    noise = Zn(1:length(dn4));
    
    % recovered signal
    yn2 = Eg*(dn4) + (noise);
    
    rx1 = zeros([1 length(yn2)]);
    for j = 1:length(yn2)
        if angle(yn2(j)) < 0
            ang = angle(yn2(j)) + 2*pi;
        else
            ang = angle(yn2(j));
        end % if
        
        if (ang >= (7*pi/4)) && (ang < (pi/4))
            rx1(j) = exp(1i*0);
        elseif (ang >= (pi/4)) && (ang < (3*pi/4))
            rx1(j) = exp(1i*pi/2);
        elseif (ang >= (3*pi/4)) && (ang < (5*pi/4))
            rx1(j) = exp(1i*pi);
        elseif (ang >= (5*pi/4)) && (ang < (7*pi/4))
            rx1(j) = exp(1i*3*pi/2);
        end % if
    end % for
    
    decodedBit = [];
    for j = 1:length(yn2)
        if rx1(j) == exp(1i*0)
            decodedBit = [decodedBit, [0,0]];
        elseif rx1(j) == exp(1i*pi/2)
            decodedBit = [decodedBit, [0,1]];
        elseif rx1(j) == exp(1i*pi)
            decodedBit = [decodedBit, [1,1]];
        elseif rx1(j) == exp(1i*3*pi/2)
            decodedBit = [decodedBit, [1,0]];
        end % if
    end % for
    
    correct = 0;
    errors = 0;
    for j = 1:length(rx1)
        if decodedBit(j) == bit_seq(j)
            correct = correct + 1;
        elseif decodedBit(j) ~= bit_seq(j)
            errors = errors + 1;
        end % if
    end % for
    
    total = errors + correct;
    BER2 = [BER2 errors/total];
    
    Pb2 = [Pb2 qfunc(sqrt(2*Eb_No))];
    EbNo2 = [EbNo2 pow2db(Eb_No)];
    
%     if SNR == 3
%         % plots
%         figure(3);
%         plot(yn2, '*');
%         title('Scatterplot of Symbols at SNR=3dB  (QPSK)');
%         xlabel('Real');
%         ylabel('Imaginary');
%     end % if
    
end % for

% plots
figure(3);
plot(yn2, '*');
title('Scatterplot of Symbols at SNR=10dB  (QPSK)');
xlabel('Real');
ylabel('Imaginary');

figure(4);
semilogy(EbNo2, BER2, EbNo2, Pb2);
title('BER vs. Eb/No  (QPSK)');
xlabel('Eb/No (dB)');
ylabel('BER (Log Scale)');
legend('Actual', 'Theoretical');

%% DBPSK

M = 2; % binary
Eb = (1/2)*(1/log2(M))*Eg; % bit energy

ddn2 = [];
ddn2(1) = dn2(1);
for i = 2:length(dn2)
    ddn2(i) = dn2(i)*ddn2(i-1);
end % for

Pb3 = [];
EbNo3 = [];
BER3 = [];
for SNR = (0:1:10)
    Eb_No = db2pow(SNR); % Eb/No
    No = Eb/Eb_No; % variance (No*Eg)
    sigma = sqrt(No);
    
    Zc = sigma * x; % real part of noise
    Zs = sigma * y; % imaginary part of noise
    Zn = Zc + 1i*Zs; % total noise
    
    % recovered signal
    yn3 = Eg*(ddn2) + (Zn);
    
    temp = [];
    for i = 2:length(dn2)
        temp(i) = yn3(i)*conj(yn3(i-1));
    end % for
    
    decodedBit = zeros([1 length(yn3)]);
    for j = 1:length(yn3)
        if temp(j) > 0
            decodedBit(j) = 0;
        elseif temp(j) < 0
            decodedBit(j) = 1;
        end % if
    end % for
    
    correct = 0;
    errors = 0;
    for j = 1:length(yn3)
        if temp(j)/dn2(j) == bit_seq(j)
            correct = correct + 1;
        elseif temp(j)/dn2(j) ~= bit_seq(j)
            errors = errors + 1;
        end % if
    end % for
    
    total = errors + correct;
    BER3 = [BER3 errors/total];
    
    Pb3 = [Pb3 qfunc(sqrt(2*Eb_No))];
    EbNo3 = [EbNo3 pow2db(Eb_No)];
    
%     if SNR == 3
%         % plots
%         figure(5);
%         plot(yn3, '*');
%         title('Scatterplot of Symbols at SNR=3dB  (DBPSK)');
%         xlabel('Real');
%         ylabel('Imaginary');
%     end % if
    
end % for

% plots
figure(5);
plot(yn3, '*');
title('Scatterplot of Symbols at SNR=10dB  (DBPSK)');
xlabel('Real');
ylabel('Imaginary');

figure(6);
semilogy(EbNo3, BER3, EbNo3, Pb3);
title('BER vs. Eb/No  (DBPSK)');
xlabel('Eb/No (dB)');
ylabel('BER (Log Scale)');
legend('Actual', 'Theoretical');

%% DQPSK

M = 4; % quad
Eb = (1/2)*(1/log2(M))*Eg; % bit energy

Pb4 = [];
EbNo4 = [];
BER4 = [];
for SNR = (0:1:10)
    Eb_No = db2pow(SNR); % Eb/No
    No = Eb/Eb_No; % variance (No*Eg)
    sigma = sqrt(No);
    
    Zc = sigma * x; % real part of noise
    Zs = sigma * y; % imaginary part of noise
    Zn = Zc + 1i*Zs; % total noise
    noise = Zn(1:length(dn4));
    
    % recovered signal
    yn4 = Eg*(dn4) + (noise);
    
    rx2 = zeros([1 length(yn4)]);
    for j = 1:length(yn4)
        if angle(yn4(j)) < 0
            ang = angle(yn4(j)) + 2*pi;
        else
            ang = angle(yn4(j));
        end % if
        
        if (ang >= (7*pi/4)) && (ang < (pi/4))
            rx2(j) = exp(1i*0);
        elseif (ang >= (pi/4)) && (ang < (3*pi/4))
            rx2(j) = exp(1i*pi/2);
        elseif (ang >= (3*pi/4)) && (ang < (5*pi/4))
            rx2(j) = exp(1i*pi);
        elseif (ang >= (5*pi/4)) && (ang < (7*pi/4))
            rx2(j) = exp(1i*3*pi/2);
        end % if
    end % for
    
    decodedBit = [];
    for j = 1:length(yn4)
        if rx2(j) == exp(1i*0)
            decodedBit = [decodedBit, [0,0]];
        elseif rx2(j) == exp(1i*pi/2)
            decodedBit = [decodedBit, [0,1]];
        elseif rx2(j) == exp(1i*pi)
            decodedBit = [decodedBit, [1,1]];
        elseif rx2(j) == exp(1i*3*pi/2)
            decodedBit = [decodedBit, [1,0]];
        end % if
    end % for
    
    correct = 0;
    errors = 0;
    for j = 1:length(rx2)
        if decodedBit(j) == bit_seq(j)
            correct = correct + 1;
        elseif decodedBit(j) ~= bit_seq(j)
            errors = errors + 1;
        end % if
    end % for
    
    total = errors + correct;
    BER4 = [BER4 errors/total];
    
    Pb4 = [Pb4 qfunc(sqrt(2*Eb_No))];
    EbNo4 = [EbNo4 pow2db(Eb_No)];
end % for

% plots
figure(7);
plot(yn4, '*');
title('Scatterplot of Symbols at SNR=10dB  (DQPSK)');
xlabel('Real');
ylabel('Imaginary');

figure(8);
semilogy(EbNo4, BER4, EbNo4, Pb4);
title('BER vs. Eb/No  (DQPSK)');
xlabel('Eb/No (dB)');
ylabel('BER (Log Scale)');
legend('Actual', 'Theoretical');

figure(9);
semilogy(EbNo1,BER1, EbNo2,BER2, EbNo3,BER3, EbNo4,BER4);
title('BER vs. Eb/No');
xlabel('Eb/No (dB)');
ylabel('BER (Log Scale)');
legend('BPSK', 'QPSK', 'DBPSK', 'DQPSK');

% figure(10);
% semilogy(EbNo1,Pb1, EbNo2,Pb2, EbNo3,Pb3, EbNo4,Pb4);
% title('BER vs. Eb/No');
% xlabel('Eb/No (dB)');
% ylabel('BER (Log Scale)');
% legend('BPSK', 'QPSK', 'DBPSK', 'DQPSK');
