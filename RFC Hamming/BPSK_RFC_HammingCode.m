%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          RFC & Hamming Code BPSK                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global Variables
SNR_dB = 0:5:50;
SNR = 10.^(SNR_dB/10);
Eb = 1;
std_dev = Eb./sqrt(2.*SNR);
a = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      BPSK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_BPSK = 2;
k_BPSK = log2(M_BPSK);

n = 2^a - 1;          %codeword length
k = 2^a - a - 1;      %message length

%Number of bits we want to send
numBits_BPSK = k*k_BPSK;

%Generate a message with k bits
dataIn_BPSK = randi([0 1], k , k_BPSK);

%Hamming Code
encoded_dataIn_BPSK_b0 = encode(dataIn_BPSK(:,1),n,k,'hamming/binary');

%Symbol is one bit, no need to loop
data_BPSK = encoded_dataIn_BPSK_b0;

%Modulate the Symbols
y_BPSK = pskmod(data_BPSK,M_BPSK,pi/M_BPSK,'gray');

z_BPSK = zeros(length(SNR_dB),length(y_BPSK));
%Add Noise to modulated symbols
for b = 1:length(SNR_dB)
    h = 1/sqrt(2)*(randn(1,n) + 1i*randn(1,n));
    for j = 1:length(y_BPSK)
        a = y_BPSK(j)*h(j);
        z_BPSK(b,j) = (awgn(a, SNR_dB(b)))/h(j);
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
    
    dataOut_dec = decode(dataOut(:,1),n,k,'Hamming/binary');
    
    nErrors = biterr(dataIn_BPSK,dataOut_dec);
    numErrs_BPSK(i) =  nErrors;
end

%BER Calculation
BER_BPSK = numErrs_BPSK./numBits_BPSK;
berTheory_BPSK = berfading(SNR_dB,'psk',2,1);

% BER vs. SNR Plot
figure(1);
semilogy(SNR_dB,BER_BPSK,'*');
hold on;
semilogy(SNR_dB,berTheory_BPSK);
grid
title('BER vs. SNR for BPSK');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
