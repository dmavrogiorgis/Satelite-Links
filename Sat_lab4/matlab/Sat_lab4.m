close all;
clear;
clc;

N=128;
M=2*N-1;
iter=10e3;
e=[1/5 1/6 1/8 1/10 1/20];

H = -e.*log2(e)-(1-e).*log2(1-e);
Cbsc = 1 - H;
y=zeros(M,1);
even_b=zeros(N-1,1);

bit_num=zeros(1,length(e));
BER=zeros(1,length(e));

for i=1:length(e)
    for j=1:iter
        m = randn(N,1)>0;
        b=zeros(M,1);
        for n=1:(length(m)-1)
            even_b(n) = xor(m(n),m(n+1));
        end
        b(1:2:end)=m;
        b(2:2:end)=even_b;
        x = rand(M,1);

        for k=1:length(y)
            if x(k)<e(i)
                y(k)=~b(k);
            else
                y(k)=b(k);
            end
        end
        y_Viterbi=Viterbi_alg(N,e(i),y);
        bit_num(i)=sum(m~=y_Viterbi');
        BER(i) = BER(i) + sum(m~=y_Viterbi');
    end    
    BER(i)=BER(i)/(N*iter);
end

disp('--------------------------------------------------------------------');
disp(['Probability å : ' num2str(e)]);
disp(['BER : ' num2str(BER)]);
disp('--------------------------------------------------------------------');

figure(1)
loglog(e,BER,' b-*');
title('BER of BSC');
xlabel('Error Probability å');
ylabel('BER - Bit Error Rate');
grid on;
%print -depsc lab4eps