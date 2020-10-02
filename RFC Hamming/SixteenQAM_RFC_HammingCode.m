%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        RFC & Hamming Code 16QAM                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Global Variables
SNR_dB = 0:5:50;
SNR = 10.^(SNR_dB/10);
Eb = 1;
std_dev = Eb./sqrt(2.*SNR);
a = 12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      16QAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_16QAM = 16;
k_16QAM = log2(M_16QAM);

n = 2^a - 1;          %codeword length
k = 2^a - a - 1;      %message length

%Generate a message with k bits
dataIn_16QAM = randi([0 1],k,k_16QAM);

numBits_16QAM = k_16QAM*k;

%Hamming Code
encoded_dataIn_16QAM_b0 = encode(dataIn_16QAM(:,1),n,k,'hamming/binary');
encoded_dataIn_16QAM_b1 = encode(dataIn_16QAM(:,2),n,k,'hamming/binary');
encoded_dataIn_16QAM_b2 = encode(dataIn_16QAM(:,3),n,k,'hamming/binary');
encoded_dataIn_16QAM_b3 = encode(dataIn_16QAM(:,4),n,k,'hamming/binary');

%Convert Bits to Symbols
data_16QAM = zeros(1,n);
for i = 1:n
    b0 = encoded_dataIn_16QAM_b0(i);
    b1 = encoded_dataIn_16QAM_b1(i);
    b2 = encoded_dataIn_16QAM_b2(i);
    b3 = encoded_dataIn_16QAM_b3(i);
    
    if(b0 == 0 && b1 == 0 && b2 == 0 && b3 == 0)
        data_16QAM(i) = 0;
    elseif(b0 == 0 && b1 == 0 && b2 == 0 && b3 == 1)
        data_16QAM(i) = 1;
    elseif(b0 == 0 && b1 == 0 && b2 == 1 && b3 == 0)
        data_16QAM(i) = 2;
    elseif(b0 == 0 && b1 == 0 && b2 == 1 && b3 == 1)
        data_16QAM(i) = 3;
        
    elseif(b0 == 0 && b1 == 1 && b2 == 0 && b3 == 0)
        data_16QAM(i) = 4;
    elseif(b0 == 0 && b1 == 1 && b2 == 0 && b3 == 1)
        data_16QAM(i) = 5;
    elseif(b0 == 0 && b1 == 1 && b2 == 1 && b3 == 0)
        data_16QAM(i) = 6;
    elseif(b0 == 0 && b1 == 1 && b2 == 1 && b3 == 1)
        data_16QAM(i) = 7;
        
    elseif(b0 == 1 && b1 == 0 && b2 == 0 && b3 == 0)
        data_16QAM(i) = 8;
    elseif(b0 == 1 && b1 == 0 && b2 == 0 && b3 == 1)
        data_16QAM(i) = 9;
    elseif(b0 == 1 && b1 == 0 && b2 == 1 && b3 == 0)
        data_16QAM(i) = 10;
    elseif(b0 == 1 && b1 == 0 && b2 == 1 && b3 == 1)
        data_16QAM(i) = 11;
    
    elseif(b0 == 1 && b1 == 1 && b2 == 0 && b3 == 0)
        data_16QAM(i) = 12;
    elseif(b0 == 1 && b1 == 1 && b2 == 0 && b3 == 1)
        data_16QAM(i) = 13;
    elseif(b0 == 1 && b1 == 1 && b2 == 1 && b3 == 0)
        data_16QAM(i) = 14;
    else
        data_16QAM(i) = 15;
    end
end

data_16QAM = transpose(data_16QAM);

%Modulate the Symbols
y_16QAM = qammod(data_16QAM,M_16QAM);

%Add Noise to modulated symbols
z_16QAM = zeros(length(SNR_dB),length(y_16QAM));
for b = 1:length(SNR_dB)
    h = 1/sqrt(2)*(randn(1,n) + 1i*randn(1,n));
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
    
    dataOut_dec = decode(dataOut(:,4),n,k,'Hamming/binary');
    dataOut_dec(:,2) = decode(dataOut(:,3),n,k,'Hamming/binary');
    dataOut_dec(:,3) = decode(dataOut(:,2),n,k,'Hamming/binary');
    dataOut_dec(:,4) = decode(dataOut(:,1),n,k,'Hamming/binary');
    
    nErrors = biterr(dataIn_16QAM,dataOut_dec);
    numErrs_16QAM(i) =  nErrors;
end

%BER Calculation
BER_16QAM = numErrs_16QAM./numBits_16QAM;
berTheory_16QAM = berfading(SNR_dB,'qam',16,1);

% BER vs. SNR Plot
figure(1);
semilogy(SNR_dB,BER_16QAM,'*');
hold on;
semilogy(SNR_dB,berTheory_16QAM);
grid
title('BER vs. SNR for 16QAM');
legend('Estimated BER','Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
 
