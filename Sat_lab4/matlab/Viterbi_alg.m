function [y_Viterbi] = Viterbi_alg(N,e,y)

w_c=log((1-e)/e);
w_00=zeros(1,N);
w_01=zeros(1,N);
w_10=zeros(1,N);
w_11=zeros(1,N);

track_00=zeros(1,N);
track_01=zeros(1,N);
track_10=zeros(1,N);
track_11=zeros(1,N);

for n=1:1:N
    if n==1
        w_00(n)=1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,0))*w_c;
        w_01(n)=1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,1))*w_c;
        w_10(n)=1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,0))*w_c;
        w_11(n)=1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,1))*w_c;
    elseif (n>1 && n<N) 
        path1_00 = w_00(n-1) + 1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,0))*w_c;
        path2_00 = w_10(n-1) + 1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,0))*w_c;
        [w_00(n),track_00(n)] = max([path1_00 0 path2_00 0]);
        
        path1_01 = w_00(n-1) + 1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,1))*w_c;
        path2_01 = w_10(n-1) + 1*(y(2*n-1)==0)*w_c + 1*(y(2*n)==xor(0,1))*w_c;
        [w_01(n),track_01(n)] = max([path1_01 0 path2_01 0]);
        
        path1_10 = w_01(n-1) + 1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,0))*w_c;
        path2_10 = w_11(n-1) + 1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,0))*w_c;
        [w_10(n),track_10(n)] = max([0 path1_10 0 path2_10]);
        
        path1_11 = w_01(n-1) + 1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,1))*w_c;
        path2_11 = w_11(n-1) + 1*(y(2*n-1)==1)*w_c + 1*(y(2*n)==xor(1,1))*w_c;
        [w_11(n),track_11(n)] = max([0 path1_11 0 path2_11]);    
    else
        path1_00 = w_00(n-1) + 1*(y(2*n-1)==0)*w_c;
        path2_00 = w_10(n-1) + 1*(y(2*n-1)==0)*w_c;
        [w_00(n),track_00(n)] = max([path1_00 0 path2_00 0]);
        
        path1_10 = w_01(n-1) + 1*(y(2*n-1)==1)*w_c;
        path2_10 = w_11(n-1) + 1*(y(2*n-1)==1)*w_c;
        [w_10(n),track_10(n)] = max([0 path1_10 0 path2_10]);
    end
end
path=zeros(1,N);
[~,path(N)] = max([w_00(n) 0 w_10(n) 0]);

y_Viterbi=zeros(1,N);
if path(n)== 1 
    y_Viterbi(n)=0;
elseif path(n)== 3
    y_Viterbi(n)=1;
end

for n=N:-1:1
    if n<N
        if path(n+1) == 1 
            path(n)=track_00(n+1);
        elseif path(n+1) == 2
            path(n)=track_01(n+1);
        elseif path(n) == 3
            path(n)=track_10(n+1);
        else
            path(n)=track_11(n+1);
        end
    end
    if path(n)==1
        y_Viterbi(n)=0;
    elseif path(n)==2
        y_Viterbi(n)=0;
    elseif path(n)==3
        y_Viterbi(n)=1;
    elseif path(n)==4
        y_Viterbi(n)=1;
    end
end