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
%                      BPSK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_BPSK = 2;
k_BPSK = log2(M_BPSK);
numBits_BPSK = N*k_BPSK;

%Symbol Generation
data_BPSK = randi([0 M_BPSK-1], N, 1);

%Convert Symbols to Bits
dataIn_BPSK = de2bi(data_BPSK);

%Modulate the Symbols
y_BPSK = pskmod(data_BPSK,M_BPSK,pi/M_BPSK,'gray');

z_BPSK = zeros(length(SNR_dB),length(y_BPSK));
%Add Noise to modulated symbols
for i = 1:length(SNR_dB);
    for j = 1:length(y_BPSK);
    z_BPSK(i,j) = awgn(y_BPSK(j),SNR_dB(i));
    end
end


%LLR Calculation
LLR_BPSK = zeros(length(SNR_dB),length(z_BPSK(i,:)));
for i = 1:length(SNR_dB)
    for j = 1:length(z_BPSK(i,:))
        yim = imag(z_BPSK(i,j));
        LLR_BPSK(i,j) = (2/std_dev(i)^2)*yim;
    end
end

%Transpose in order to demodulate
z_BPSK_T = transpose(z_BPSK);

%Demodulate, Calcuate Number of Errors and Bits
numErrs_BPSK = zeros(1,length(SNR_dB));
z_BPSK_demod = zeros(length(y_BPSK),length(SNR_dB));
for i = 1:length(SNR_dB)
    numErrs_BPSK(i) = 0;
    z_BPSK_demod(:,i) = pskdemod(z_BPSK_T(:,i),M_BPSK,pi/M_BPSK,'gray');
    dataOut = de2bi(z_BPSK_demod(:,i));
    nErrors = biterr(dataIn_BPSK,dataOut);
    numErrs_BPSK(i) =  nErrors;
end

%BER Calculation
BER_BPSK = numErrs_BPSK./numBits_BPSK;
berTheory_BPSK = berawgn(SNR_dB,'psk',M_BPSK,'nondiff');

% BER vs. SNR Plot
figure(figure_num);
semilogy(SNR_dB,BER_BPSK,'*');
hold on;
semilogy(SNR_dB,berTheory_BPSK);
grid
title('BER vs. SNR for BPSK');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');

figure_num = figure_num + 1;

%Plot of Symbols w/o Noise
figure(figure_num);
scatterplot(y_BPSK);
title('BPSK Without Noise');

figure_num = figure_num + 1;
j = 1;
snr = 0;

% Plot of Received Noisy Symbols for Different SNR values
for i = figure_num:1:(figure_num + graph_num)
    figure(i);
    scatterplot(z_BPSK(j,:));
    title(['BPSK - SNR = ' num2str(snr)]);
    j = j+1;
    snr = snr + SNR_step;
end

figure_num = figure_num + graph_num + 1;