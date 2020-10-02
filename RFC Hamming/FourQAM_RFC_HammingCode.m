%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         RFC & Hamming Code 4QAM                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global Variables
SNR_dB = 0:5:50;
SNR = 10.^(SNR_dB/10);
Eb = 1;
std_dev = Eb./sqrt(2.*SNR);
a = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      4QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_4QAM = 4;
k_4QAM = log2(M_4QAM);

n = 2^a - 1;          %codeword length
k = 2^a - a - 1;      %message length

%Generate a message with k bits
dataIn_4QAM = randi([0 1],k,k_4QAM);

numBits_4QAM = k_4QAM*k;

%Hamming Code
encoded_dataIn_4QAM_b0 = encode(dataIn_4QAM(:,1),n,k,'hamming/binary');
encoded_dataIn_4QAM_b1 = encode(dataIn_4QAM(:,2),n,k,'hamming/binary');

%Convert Bits to Symbols
data_4QAM = zeros(1,n);
for i = 1:n
    b0 = encoded_dataIn_4QAM_b0(i);
    b1 = encoded_dataIn_4QAM_b1(i);
    
    if (b0 == 0 && b1 == 0)
        data_4QAM(i) = 0;
    elseif(b0 == 0 && b1 == 1)
        data_4QAM(i) = 1;
    elseif(b0 == 1 && b1 == 0)
        data_4QAM(i) = 2;
    elseif(b0 == 1 && b1 == 1)
        data_4QAM(i) = 3;
    end
end

data_4QAM = transpose(data_4QAM);

%Modulate the Symbols
y_4QAM = qammod(data_4QAM,M_4QAM);

%Add Noise to the modulated symbols
z_4QAM = zeros(length(SNR_dB),length(y_4QAM));
for b = 1:length(SNR_dB)
    h = 1/sqrt(2)*(randn(1,n) + 1i*randn(1,n));
    for j = 1:length(y_4QAM)
        a = y_4QAM(j)*h(j);
        z_4QAM(b,j) = (awgn(a ,SNR_dB(b)))/h(j);
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
    
    dataOut_dec = decode(dataOut(:,2),n,k,'Hamming/binary');
    dataOut_dec(:,2) = decode(dataOut(:,1),n,k,'Hamming/binary');
    
    nErrors = biterr(dataIn_4QAM,dataOut_dec);
    numErrs_4QAM(i) = nErrors;
end


%BER Calculation
BER_4QAM = numErrs_4QAM./numBits_4QAM;
berTheory_4QAM = berfading(SNR_dB,'qam',4,1);

% BER vs. SNR Plot
figure(1);
semilogy(SNR_dB,BER_4QAM,'*');
hold on;
semilogy(SNR_dB,berTheory_4QAM);
grid
title('BER vs. SNR for 4QAM');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
 

