function J=tdoaFTdpd(Rk,B,P_bs,P_ue2,D,fn,lamda)
c=3e8;
K=size(Rk,3);
Nr=size(P_bs,1);
M=size(B,1);
Kvec=-2*pi/lamda*(([P_ue2,4]-P_bs).'./sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).')).').';
A=zeros(M,Nr);
for j=1:Nr
A(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau=sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).'))/c;%zeros(1,Nr);

Dk=zeros(Nr);
for k=1:K
    Bk=[];
    for j=1:Nr
        Bk=blkdiag(Bk,B*A(:,j)*exp(-2j*pi*fn(k)*tau(j)));
    end
    Dk=Dk+Bk'*Rk(:,:,k)*Bk;
end
    [~,te,~]=svd(Dk);
    J=-te(1);
end