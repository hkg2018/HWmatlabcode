function J=MVDR_dpd(invRy,Xph,B,P_bs,P_ue2,D,fnph,lamda)
c=3e8;
N_bs=size(P_bs,1);
M=size(B,1);
Kvec=-2*pi/lamda*(([P_ue2,4]-P_bs).'./sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).')).').';
a=zeros(M,N_bs);
for nbs=1:N_bs
    a(:,nbs)=exp(1j*(Kvec(nbs,:)*D(:,:,nbs)).');
end
tau0=sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).'))/c;%zeros(1,Nr);
J=0;
for nbs=1:N_bs
    
    Aj=kron(diag(exp(-2j*pi*fnph*tau0(nbs)))*Xph,B*a(:,nbs));
    J=J+Aj'*invRy(:,:,nbs)*Aj;
end
J=-1/J;
end
