%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          RFC Code - 16QAM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%                      16QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_16QAM = 16;
k_16QAM = log2(M_16QAM);
numBits_16QAM = N*k_16QAM;

%Symbol Generation
data_16QAM = randi([0 M_16QAM-1], N, 1);

%Convert Symbols to Bits
dataIn_16QAM = de2bi(data_16QAM);

%Modulate the Symbols
y_16QAM = qammod(data_16QAM,M_16QAM);

%Add Noise to modulated symbols
z_16QAM = zeros(length(SNR_dB),length(y_16QAM));
for b = 1:length(SNR_dB)
    h = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    for j = 1:length(y_16QAM)
        a = y_16QAM(j)*h(j);
        z_16QAM(b,j) = (awgn(a, SNR_dB(b)))/h(j);
    end
end

%Transpose in order to demodulate
z_16QAM_T = transpose(z_16QAM);

%Demodulate, Calcuate Number of Errors and Bits
numErrs_16QAM = zeros(1,length(SNR_dB));
z_16QAM_demod = zeros(length(y_16QAM),length(SNR_dB));
for i = 1:length(SNR_dB)
    numErrs_16QAM(i) = 0;
    z_16QAM_demod(:,i) = qamdemod(z_16QAM_T(:,i),M_16QAM);
    dataOut = de2bi(z_16QAM_demod(:,i));
    nErrors = biterr(dataIn_16QAM,dataOut);
    numErrs_16QAM(i) =  nErrors;
end

%BER Calculation
BER_16QAM = numErrs_16QAM./numBits_16QAM;
berTheory_16QAM = berfading(SNR_dB,'qam',16,1);

% BER vs. SNR Plot
figure(figure_num);
semilogy(SNR_dB,BER_16QAM,'*');
hold on;
semilogy(SNR_dB,berTheory_16QAM);
grid
title('BER vs. SNR for 16QAM');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
 
figure_num = figure_num + 1;

%Plot of Symbols w/o Noise
figure(figure_num);
scatterplot(y_16QAM);
title('16QAM Without Noise');

figure_num = figure_num + 1;
j = 1;
snr = 0;

% Plot of Received Noisy Symbols for Different SNR values
for i = figure_num:1:(figure_num + graph_num)
    figure(i);
    scatterplot(z_16QAM(j,:));
    title(['16QAM - SNR = ' num2str(snr)]);
    j = j+1;
    snr = snr + SNR_step;
    xlim([-4 4]);
    ylim([-4 4]);
end

figure_num = figure_num + graph_num + 1;