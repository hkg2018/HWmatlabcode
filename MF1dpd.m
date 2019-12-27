function J=MF1dpd(Ydpd_t,B,P_bs,P_ue_xy,P_ue_z,D,fn0,lamda,h_re,tau_re)
c=3e8;
Nr=size(P_bs,1);
M=size(B,1);
Kvec=-2*pi/lamda*(([P_ue_xy,4]-P_bs).'./sqrt(diag(([P_ue_xy,4]-P_bs)*([P_ue_xy,4]-P_bs).')).').';
a=zeros(M,Nr);
for j=1:Nr
    a(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau0=sqrt(diag(([P_ue_xy,4]-P_bs)*([P_ue_xy,4]-P_bs).'))/c;%zeros(1,Nr);
J=0;
for j=1:Nr
    
    Aj=kron(exp(-2j*pi*fn0*tau0(j)),B*a(:,j));
    
    % J=J-abs(sum(diag(Y(:,:,j)'*Aj)))^2;
    J=J+1/abs(sum(diag(Ydpd_t(:,:,j)'*Aj)))^2;
end
J=-1/J;
end
