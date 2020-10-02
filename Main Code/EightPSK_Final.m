%Global Variables
SNR_dB = 0:5:50;
SNR = 10.^(SNR_dB/10);
Eb = 1;
std_dev = Eb./sqrt(2.*SNR);
N = 100000;
figure_num = 1;
graph_num = 10;
SNR_step = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      8PSK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_8PSK = 8;
k_8PSK = log2(M_8PSK);
numBits_8PSK = N*k_8PSK;

%Symbol Generation
data_8PSK = randi([0 M_8PSK-1], N, 1);

%Convert Symbols to Bits
dataIn_8PSK = de2bi(data_8PSK);

%Modulate the Symbols
y_8PSK = pskmod(data_8PSK,M_8PSK,pi/M_8PSK,'gray');

z_8PSK = zeros(length(SNR_dB),length(y_8PSK));
%Add Noise to modulated symbols
for i = 1:length(SNR_dB);
    for j = 1:length(y_8PSK);
    z_8PSK(i,j) = awgn(y_8PSK(j),SNR_dB(i));
    end
end

%Transpose in order to demodulate
z_8PSK_T = transpose(z_8PSK);

%Demodulate, Calcuate Number of Errors and Bits
numErrs_8PSK = zeros(1,length(SNR_dB));
z_8PSK_demod = zeros(length(y_8PSK),length(SNR_dB));
for i = 1:length(SNR_dB)
    numErrs_8PSK(i) = 0;
    z_8PSK_demod(:,i) = pskdemod(z_8PSK_T(:,i),M_8PSK,pi/M_8PSK,'gray');
    dataOut = de2bi(z_8PSK_demod(:,i));
    nErrors = biterr(dataIn_8PSK,dataOut);
    numErrs_8PSK(i) = nErrors;
end

%BER Calculation
BER_8PSK = numErrs_8PSK./numBits_8PSK;
berTheory_8PSK = berawgn(SNR_dB,'psk',M_8PSK,'nondiff');

% BER vs. SNR Plot
figure(figure_num);
semilogy(SNR_dB,BER_8PSK,'*');
hold on;
semilogy(SNR_dB,berTheory_8PSK);
grid
title('BER vs. SNR for 8PSK');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');

figure_num = figure_num + 1;

%Plot of Symbols w/o Noise
figure(figure_num);
scatterplot(y_8PSK);
title('8PSK Without Noise');

figure_num = figure_num + 1;
j = 1;
snr = 0;

% Plot of Received Noisy Symbols for Different SNR values
for i = figure_num:1:(figure_num + graph_num)
    figure(i);
    scatterplot(z_8PSK(j,:));
    title(['8PSK - SNR = ' num2str(snr)]);
    j = j+1;
    snr = snr + SNR_step;
end

figure_num = figure_num + graph_num + 1;