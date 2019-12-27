function J=MF0dpd(Ydpd_t,B,P_bs,P_ue_xy,P_ue_z,D,fn,lamda,h_re,tau_re)
c=3e8;
P_ue=[P_ue_xy,P_ue_z];
N_bs=size(P_bs,1);
N=length(fn);
M=size(B,2);
Kvec=-2*pi/lamda*((P_ue-P_bs).'./sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).')).').';
a=zeros(M,N_bs);
for j=1:N_bs
    a(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau0=sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).'))/c;%zeros(1,Nr);
Phi=zeros(N);
for j=1:Nr
    
    Aj=diag(h_re(:,j).*exp(2j*pi*fn*(tau_re(j)-tau0(j))));%kron(eye(N),B*a(:,j))*
    
    Phi=Phi+Aj'*Ydpd_t(:,j)*Ydpd_t(:,j)'*Aj;
end
[~,qyz,~]=svd(Phi);


J=-max(abs(diag(qyz)))^2;
end
