close all;
clear all;
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
disp('1st Question');
disp('------------');
disp(['SNR(dB) : ' num2str(SNR_dB(1))]);
disp(['Experimental BER : ' num2str(BER(1))]);
disp(['Theoritical BER : ' num2str(BER_theor(1))]);
disp(['Accurate Theoritical BER : ' num2str(BER_theor1(1))]);
disp('-----------------------------------');

disp('2nd Question');
disp('------------');
disp(['SNR(dB) : ' num2str(SNR_dB(2:end))]);
disp(['Theoritical BER : ' num2str(BER(2:end))]);
disp(['Theoritical BER : ' num2str(BER_theor(2:end))]);
disp(['Accurate Theoritical BER : ' num2str(BER_theor1(2:end))]);
disp('-----------------------------------------------------------------------------------------------------------');

figure(1)
semilogy(SNR_dB,BER,'-*') 
hold on;
semilogy(SNR_dB,BER_theor,'-o')
semilogy(SNR_dB,BER_theor1,'-X')
xlabel('SNR (dB)')
ylabel('BER - Bit Error Rate')
title('MSK BER through OQPSK')
legend('Experimental BER', 'Theoretical BER','Theoretical BER(More Accurate)');
grid on;
%print -depsc epsFig