function [J_Musicdpd,J_MF0dpd]=indoor_dpd(Ydpd_t,Us,B,P_bs,P_ue_xy,P_ue_z,D1,fn,lamda,h_re,tau_re,LenFreWindow)
c=3e8;
P_ue=[P_ue_xy,P_ue_z];
N_bs=size(P_bs,1);
N=length(fn);
M=size(B,2);
fnph=fn(1:LenFreWindow);
h_reph=h_re(1:LenFreWindow,:);
Kvec=-2*pi/lamda*((P_ue-P_bs).'./sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).')).').';
a=zeros(M,N_bs);
for j=1:N_bs
    a(:,j)=exp(1j*(Kvec(j,:)*D1(:,:,j)).');
end
tau0=sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).'))/c;%zeros(1,Nr);
Phi_music=zeros(LenFreWindow);
Phi_MF=zeros(N);
for j=1:N_bs
    
    Aj_music=diag(h_reph(:,j).*exp(2j*pi*fnph*(tau_re(j)-tau0(j))));%kron(eye(LenFreWindow),B*a(:,j))*
    Aj_MF=diag(h_re(:,j).*exp(2j*pi*fn*(tau_re(j)-tau0(j))));%kron(eye(N),B*a(:,j))*
    % J=J-abs(sum(diag(Y(:,:,j)'*Aj)))^2;
    Phi_music=Phi_music+Aj_music'*Aj_music-Aj_music'*(Us(:,1:N_paths(j),j)*Us(:,1:N_paths(j),j)')*Aj_music;
    Phi_MF=Phi_MF+Aj_MF'*Ydpd_t(:,j)*Ydpd_t(:,j)'*Aj_MF;
end
[~,qyz_music,~]=svd(Phi_music);
[~,qyz_MF,~]=svd(Phi_MF);

J_Musicdpd=1/min(abs(diag(qyz_music)))^2;
J_MF0dpd=max(abs(diag(qyz_MF)))^2;
end