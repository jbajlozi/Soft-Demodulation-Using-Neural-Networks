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
for i = 1:length(SNR_dB)
    for j = 1:length(y_16QAM)
        z_16QAM(i,j) = awgn(y_16QAM(j),SNR_dB(i));
    end
end

%LLR Calculation
LLR_b0_16QAM = zeros(length(SNR_dB),length(z_16QAM(i,:)));
LLR_b1_16QAM = zeros(length(SNR_dB),length(z_16QAM(i,:)));
LLR_b2_16QAM = zeros(length(SNR_dB),length(z_16QAM(i,:)));
LLR_b3_16QAM = zeros(length(SNR_dB),length(z_16QAM(i,:)));
for i = 1:length(SNR_dB)
    for j = 1:length(z_16QAM(i,:))
        yre = real(z_16QAM(i,j));
        yim = imag(z_16QAM(i,j));
        %b0
        if(yre < -2)
            LLR_b0_16QAM(i,j) = (1/std_dev(i)^2)*4*(yre+1);
        elseif(yre >= -2 && yre < 2)
            LLR_b0_16QAM(i,j) = (2/std_dev(i)^2)*yre;
        else 
            LLR_b0_16QAM(i,j) = (4/std_dev(i)^2)*(yre-1);
        end
        %b1
        if(yre < -0)
            LLR_b1_16QAM(i,j) = (2/std_dev(i)^2)*(yre+2);
        else 
            LLR_b1_16QAM(i,j) = (2/std_dev(i)^2)*(-yre+2);
        end
        
        %b2
        if(yim < -2)
            LLR_b2_16QAM(i,j) = (1/std_dev(i)^2)*4*(yim+1);
        elseif(yim >= -2 && yim < 2)
            LLR_b2_16QAM(i,j) = (2/std_dev(i)^2)*yim;
        else 
            LLR_b2_16QAM(i,j) = (4/std_dev(i)^2)*(yim-1);
        end
        %b3
        if(yim < -0)
            LLR_b3_16QAM(i,j) = (2/std_dev(i)^2)*(yim+2);
        else 
            LLR_b3_16QAM(i,j) = (2/std_dev(i)^2)*(-yim+2);
        end
        
 
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
berTheory_16QAM = berawgn(SNR_dB,'qam',M_16QAM);

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
end

figure_num = figure_num + graph_num + 1;