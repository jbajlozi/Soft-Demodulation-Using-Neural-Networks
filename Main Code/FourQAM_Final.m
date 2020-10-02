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
%                      4QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_4QAM = 4;
k_4QAM = log2(M_4QAM);
numBits_4QAM = N*k_4QAM;

%Symbol Generation
data_4QAM = randi([0 M_4QAM-1], N, 1);

%Convert Symbols to Bits
dataIn_4QAM = de2bi(data_4QAM);

%Modulate the Symbols
y_4QAM = qammod(data_4QAM,M_4QAM);

%Add Noise to the modulated symbols
z_4QAM = zeros(length(SNR_dB),length(y_4QAM));
for i = 1:length(SNR_dB)
    for j = 1:length(y_4QAM)
        z_4QAM(i,j) = awgn(y_4QAM(j),SNR_dB(i));
    end
end


%LLR Calculation
LLR_BPSK_b0 = zeros(length(SNR_dB),length(z_4QAM(1,:)));
LLR_BPSK_b1 = zeros(length(SNR_dB),length(z_4QAM(1,:)));
for i = 1:length(SNR_dB)
    for j = 1:length(z_4QAM(i,:))
        yre = real(z_4QAM(i,j));
        yim = imag(z_4QAM(i,j));
        %b0
        LLR_BPSK_b0(i,j) = (2/std_dev(i)^2)*yre;
        %b1
        LLR_BPSK_b1(i,j) = -(2/std_dev(i)^2)*yim;
    end
end

%Transpose in order to demodulate
z_4QAM_T = transpose(z_4QAM);

%Demodulate, Calcuate Number of Errors and Bits
z_4QAM_demod = zeros(length(y_4QAM),length(SNR_dB));
numErrs_4QAM = zeros(1,length(SNR_dB));
for i = 1:length(SNR_dB)
    numErrs_4QAM(i) = 0;
    z_4QAM_demod(:,i) = qamdemod(z_4QAM_T(:,i),M_4QAM);
    dataOut = de2bi(z_4QAM_demod(:,i));
    nErrors = biterr(dataIn_4QAM,dataOut);
    numErrs_4QAM(i) = nErrors;
end


%BER Calculation
BER_4QAM = numErrs_4QAM./numBits_4QAM;
berTheory_4QAM = berawgn(SNR_dB,'qam',M_4QAM);

% BER vs. SNR Plot
figure(figure_num);
semilogy(SNR_dB,BER_4QAM,'*');
hold on;
semilogy(SNR_dB,berTheory_4QAM);
grid
title('BER vs. SNR for 4QAM');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
 
figure_num = figure_num + 1;

%Plot of Symbols w/o Noise
figure(figure_num);
scatterplot(y_4QAM);
title('4QAM - Without Noise');

figure_num = figure_num + 1;
j = 1;
snr = 0;

% Plot of Received Noisy Symbols for Different SNR values
for i = figure_num:1:(figure_num + graph_num)
    figure(i);
    scatterplot(z_4QAM(j,:));
    title(['4QAM - SNR = ' num2str(snr)]);
    j = j+1;
    snr = snr + SNR_step;
end
