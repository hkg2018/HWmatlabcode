function J=Music_dpd(Us,N_paths,B,P_bs,P_ue_xy,P_ue_z,D,fnph,lamda,h_re,tau_re)
c=3e8;
P_ue=[P_ue_xy,P_ue_z];
N_bs=size(P_bs,1);
LenFreWindow=length(fnph);
M=size(B,2);
Kvec=-2*pi/lamda*((P_ue-P_bs).'./sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).')).').';
a=zeros(M,N_bs);
for j=1:N_bs
    a(:,j)=exp(1j*(Kvec(j,:)*D(:,:,j)).');
end
tau0=sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).'))/c;%zeros(1,Nr);
Phi=zeros(LenFreWindow);
for j=1:N_bs
    
    Aj=diag(h_re(1:LenFreWindow,j).*exp(2j*pi*fnph*(tau_re(j)-tau0(j))));%kron(eye(N),B*a(:,j))*
    
    % J=J-abs(sum(diag(Y(:,:,j)'*Aj)))^2;
    Phi=Phi+Aj'*Aj-Aj'*(Us(:,1:N_paths(j),j)*Us(:,1:N_paths(j),j)')*Aj;
end
[~,qyz,~]=svd(Phi);

J=-1/min(abs(diag(qyz)));
end
