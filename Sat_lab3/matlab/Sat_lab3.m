close all;
clear;
clc

N=10^5;
A=1;
T=10^-3;
SNR_dB=5:12;
SNR=10.^(SNR_dB/10);
b=A^2*T./SNR;

BER=zeros(1,length(SNR));
BER_theor=zeros(1,length(SNR));
BER_theor1=zeros(1,length(SNR));
BER_Viterbi=zeros(1,length(SNR));

out_symbols=zeros(N,1);
out_symbols_odd=zeros(N/2,1);
out_symbols_even=zeros(N/2,1);

in_bits = randn(N,1)>0; %Creation of input bits {0,1}
in_symbols = 2*in_bits-1; % [-1 1 -1 -1 1 -1 -1 -1 ].';Symbols we transmit {-1,+1}

xQ1=-1;
x1=1;

xI_n=zeros(N/2,1);
xQ_n=zeros(N/2,1);
for k=1:1:length(xI_n)
    if k==1     %Symbol x0
        xI_n(k)= -xQ1*x1;
        xQ_n(k)= -xI_n(k)*in_symbols(k);
    else    %Other symbols
        xI_n(k)= -xQ_n(k-1)*in_symbols(2*k-2);
        xQ_n(k)= -xI_n(k)*in_symbols(2*k-1);
    end
end

z_n=xI_n + 1j*xQ_n; %Construction of O-QPSK sympols

for n=1:length(SNR)
    n_n = randn(N/2,1) + 1j*randn(N/2,1); 
    y_n=A*T*z_n + sqrt(b(n)*T)*n_n;
    yI_n = sign(real(y_n)); 
    yQ_n = sign(imag(y_n));
    for m=1:1:(length(y_n)-1)
        if m==1
            out_symbols_odd(m) =  -yQ_n(m)*yI_n(m);
            out_symbols_even(m) = -yQ_n(m)*yI_n(m+1);
        else
            out_symbols_odd(m) =  -yQ_n(m)*yI_n(m);
            out_symbols_even(m) = -yQ_n(m)*yI_n(m+1);
        end 
    end
    out_symbols_odd(N/2) =  -yQ_n(N/2)*yI_n(N/2); 
    out_symbols_even(N/2) =  -yQ_n(N/2-1)*yI_n(N/2); 
    out_symbols(1:2:end)=out_symbols_odd;
    out_symbols(2:2:end)=out_symbols_even;
    
    BER(n) = sum(in_symbols~=out_symbols)/N;
    BER_theor(n) = qfunc(sqrt(SNR(n))); %B-PSK theoretical BER
    BER_theor1(n) = 2*qfunc(sqrt(SNR(n))) - (qfunc(sqrt(SNR(n))))^2; % Accyrate theoretical BER
end

s_1 = [-(2*A*sqrt(T)*1i)/pi; (A*sqrt(T)*sqrt(pi^2 - 4))/pi];
s1  = [A*sqrt(T); 0];

phi=zeros(1,N+1);
r_n=zeros(2,N);
out_symbols_Viterbi=zeros(1,N);

for n=1:length(SNR)
    n_n1 = sqrt(b(n))*(randn(1,N) + 1j*randn(1,N));
    n_n2 = sqrt(b(n))*(randn(1,N) + 1j*randn(1,N));
    phi(1)=0;
    for m=1:N
        phi(m+1)= phi(m) + in_symbols(m)*(pi/2);
        if in_symbols(m)==1
            r_n(:,m) = s1.*exp(1j*phi(m)) + [n_n1(m); n_n2(m)];
        else
            r_n(:,m) = s_1.*exp(1j*phi(m)) + [n_n1(m); n_n2(m)];
        end
    end
    out_symbols_Viterbi=Viterbi_alg(N,s1,s_1,r_n);
    BER_Viterbi(n) = sum(in_symbols~=out_symbols_Viterbi.')/N;
end

disp('1st Question');
disp('------------');
disp(['SNR(dB) : ' num2str(SNR_dB(1))]);
disp(['Experimental BER : ' num2str(BER(1))]);
disp(['Theoritical BER : ' num2str(BER_theor(1))]);
disp(['Accurate Theoritical BER : ' num2str(BER_theor1(1))]);
disp(['Viterbi Algorithm BER : ' num2str(BER_Viterbi(1))]);
disp('-----------------------------------');

disp('2nd Question');
disp('------------');
disp(['SNR(dB) : ' num2str(SNR_dB(2:end))]);
disp(['Theoritical BER : ' num2str(BER(2:end))]);
disp(['Theoritical BER : ' num2str(BER_theor(2:end))]);
disp(['Accurate Theoritical BER : ' num2str(BER_theor1(2:end))]);
disp(['Viterbi Algorithm BER : ' num2str(BER_Viterbi(2:end))]);
disp('-----------------------------------------------------------------------------------------------------------');

figure(1)
semilogy(SNR_dB,BER,'-*') 
hold on;
semilogy(SNR_dB,BER_theor,'-o')
semilogy(SNR_dB,BER_Viterbi,'-V')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate')
title('MSK BER through OQPSK')
legend('Experimental BER', 'Theoretical BER','Viterbi Aglorithm BER');
hold off;
grid on;
% print -depsc epsFig1

figure(2)
semilogy(SNR_dB,BER_Viterbi,'-V')
hold on;
semilogy(SNR_dB,BER_theor,'-o')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate')
title('MSK BER-Viterbi Algorithm')
legend('Viterbi Aglorithm BER','Theoretical BER');
grid on;
% print -depsc epsFig2

figure(3)
semilogy(SNR_dB,BER,'-*') 
hold on;
semilogy(SNR_dB,BER_theor,'-o')
semilogy(SNR_dB,BER_theor1,'-X')
semilogy(SNR_dB,BER_Viterbi,'-V')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate')
title('MSK BER-All methods')
legend('Experimental BER', 'Theoretical BER','Theoretical BER(More Accurate)','Viterbi Aglorithm BER');
hold off;
grid on;
% print -depsc epsFig3