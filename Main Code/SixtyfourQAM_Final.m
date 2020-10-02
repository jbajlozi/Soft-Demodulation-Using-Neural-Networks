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
%                      64QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_64QAM = 64;
k_64QAM = log2(M_64QAM);
numBits_64QAM = N*k_64QAM;

%Symbol Generation
data_64QAM = randi([0 M_64QAM-1], N, 1);

%Convert Symbols to Bits
dataIn_64QAM = de2bi(data_64QAM);

%Modulate the Symbols
y_64QAM = qammod(data_64QAM,M_64QAM);

%Add Noise to modulated symbols
z_64QAM = zeros(length(SNR_dB),length(y_64QAM));
for i = 1:length(SNR_dB)
    for j = 1:length(y_64QAM)
    z_64QAM(i,j) = awgn(y_64QAM(j),SNR_dB(i));
    end
end

%Transpose in order to demodulate
z_64QAM_T = transpose(z_64QAM);

%Demodulate, Calcuate Number of Errors and Bits
z_64QAM_demod = zeros(length(y_64QAM),length(SNR_dB));
numErrs_64QAM = zeros(1,length(SNR_dB));
for i = 1:length(SNR_dB)
    numErrs_64QAM(i) = 0;
    z_64QAM_demod(:,i) = qamdemod(z_64QAM_T(:,i),M_64QAM);
    dataOut = de2bi(z_64QAM_demod(:,i));
    nErrors = biterr(dataIn_64QAM,dataOut);
    numErrs_64QAM(i) = nErrors;
end

%BER Calculation
BER_64QAM = numErrs_64QAM./numBits_64QAM;
berTheory_64QAM = berawgn(SNR_dB,'qam',M_64QAM);

% BER vs. SNR Plot
figure(figure_num);
semilogy(SNR_dB,BER_64QAM,'*');
hold on;
semilogy(SNR_dB,berTheory_64QAM);
grid
title('BER vs. SNR for 64QAM');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');

figure_num = figure_num + 1;

%Plot of Symbols w/o Noise
figure(figure_num);
scatterplot(y_64QAM);
title('64QAM Without Noise');

figure_num = figure_num + 1;
j = 1;
snr = 0;

% Plot of Received Noisy Symbols for Different SNR values
for i = figure_num:1:(figure_num + graph_num)
    figure(i);
    scatterplot(z_64QAM(j,:));
    title(['64QAM - SNR = ' num2str(snr)]);
    j = j+1;
    snr = snr + SNR_step;
end

figure_num = figure_num + graph_num + 1;