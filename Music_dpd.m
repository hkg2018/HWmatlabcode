function J=Music_dpd(Us,Xph,B,P_bs,P_ue2,D,fnph,lamda)
c=3e8;
N_bs=size(P_bs,1);
M=size(B,1);
Kvec=-2*pi/lamda*(([P_ue2,4]-P_bs).'./sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).')).').';
a=zeros(M,N_bs);
for j=1:N_bs
    a(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau0=sqrt(diag(([P_ue2,4]-P_bs)*([P_ue2,4]-P_bs).'))/c;%zeros(1,Nr);
J=0;
for j=1:N_bs
    
    Aj=kron(diag(exp(-2j*pi*fnph*tau0(j)))*Xph,B*a(:,j));
    
    % J=J-abs(sum(diag(Y(:,:,j)'*Aj)))^2;
    J=J+norm(Aj)^2-norm(Aj'*Us(:,:,j))^2;
end
J=-1/J;
end
