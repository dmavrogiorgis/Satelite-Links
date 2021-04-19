function [out_symbols_Viterbi] = Viterbi_alg(N, s1, s_1, r_n)
w_0=zeros(1,N);
w_pi2=zeros(1,N);
w_pi=zeros(1,N);
w_3pi2=zeros(1,N);

track_0=zeros(1,N);
track_pi2=zeros(1,N);
track_pi=zeros(1,N);
track_3pi2=zeros(1,N);

for n=1:1:N
    if n==1
        w_3pi2(n) = real((r_n(:,1)')*s_1*exp(1i*0)); % This is for the 1st step.  
        w_pi2(n) = real((r_n(:,1)')*s1*exp(1i*0));
    elseif n~=1
        if mod(n,2)==0 % Even bits
            path_w1 = real((r_n(:,n)')*s_1*exp(1i*3*pi/2)) + w_3pi2(n-1); %from 3pi/2 to pi
            path_w2 = real((r_n(:,n)')*s1*exp(1i*pi/2)) + w_pi2(n-1); %from pi/2 to pi
            [w_pi(n),track_pi(n)] = max([path_w1 0 path_w2 0]);

            path_w3 = real((r_n(:,n)')*s1*exp(1i*3*pi/2)) + w_3pi2(n-1); %from 3pi/2 to 0
            path_w4 = real((r_n(:,n)')*s_1*exp(1i*pi/2)) + w_pi2(n-1); %from pi/2 to 0
            [w_0(n),track_0(n)] = max([path_w3 0 path_w4 0]);
            
        elseif mod(n,2)~=0 % Odd bits
            path_w1 = real((r_n(:,n)')*s1*exp(1j*pi)) + w_pi(n-1);%from pi to 3pi/2
            path_w2 = real((r_n(:,n)')*s_1) + w_0(n-1);%from 0 to 3pi/2
            [w_3pi2(n),track_3pi2(n)] = max([0 path_w1 0 path_w2]);

            path_w3 = real((r_n(:,n)')*s_1*exp(1j*pi)) + w_pi(n-1); %from pi to pi/2
            path_w4 = real((r_n(:,n)')*s1*exp(1j*0)) + w_0(n-1); %from 0 to 3pi/2
            [w_pi2(n),track_pi2(n)] = max([0 path_w3 0 path_w4]);
        end
    end
end
path=zeros(1,N+1);
path(1) = 4;
if mod(n,2)~=0 %Odd num of symbols
    [~,path(N+1)] = max([w_3pi2(n) 0 w_pi2(n) 0]);
elseif mod(n,2)==0 %Even num of symbols
    [~,path(N+1)] = max([0 w_pi(n) 0 w_0(n)]);
end

for n=N:-1:1
    if n>1
        if mod(n,2)==0 % Even symbols
            [~,index] = max([0 w_pi(n) 0 w_0(n)]);
            track = [0 track_pi(n) 0 track_0(n)];
            path(n) = track(index);
        elseif mod(n,2)~=0 % Odd symbols
            [~,index] = max([w_3pi2(n) 0 w_pi2(n) 0]);
            track = [track_3pi2(n) 0 track_pi2(n) 0];
            path(n) = track(index);
        end
    end
    if path(n)-path(n+1)==-1
        out_symbols_Viterbi(n) = -1;
    elseif path(n)-path(n+1)==1
        out_symbols_Viterbi(n) = 1;
    elseif path(n)-path(n+1)==-3
        out_symbols_Viterbi(n) = 1;
    elseif path(n)-path(n+1)==3
        out_symbols_Viterbi(n) = -1;
    end 
end
end