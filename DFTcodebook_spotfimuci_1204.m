function T=DFTcodebook_spotfimuci_1204_2(M,N)
T0=zeros(M*N,M*N);
for Numn=0:(N-1)
for Numm=0:(M-1)
       for n=0:(N-1)
        for m=0:(M-1)
            WM=m/M;
            WN=n/N;
%             if(WM>0.5)
%                 WM=WM-1;
%             end
%             if(WN>0.5)
%                 WN=WN-1;

             T0(n*M+m+1,Numn*M+Numm+1)=exp(-1i*2*pi*((Numm)*WM+(Numn)*WN));%此时需要共轭转置                         
        end
end
%% &&&&&&&&&&&&&&&&&&&&&&
    T1=T0(:,Numn*M+Numm+1)';
    if((Numn*M+Numm+1)==1)
        T3=T1;
    else
        T3=[T3;T1];
    end
end
end
T=T3;%T为M*N
end