function J=MFdpd(Y,B,P_bs,P_ue2,D,fn0,lamda)
c=3e8;
Nr=size(P_bs,1);
M=size(B,1);
Kvec=-2*pi/lamda*(([P_ue2,4]-P_bs).'./sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).')).').';
a=zeros(M,Nr);
for j=1:Nr
    a(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau0=sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).'))/c;%zeros(1,Nr);
J=0;
for j=1:Nr
    
    Aj=kron(exp(-2j*pi*fn0*tau0(j)),B*a(:,j));
    
    J=J-abs(sum(diag(Y(:,:,j)'*Aj)))^2;
end
end
