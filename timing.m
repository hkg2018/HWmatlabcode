disp('Smusic-tau')
Rtau_e=(0.1:.1:20)*1e-9;
Y_t0=permute(Y_t,[1,3,2,4]);
Us1=zeros(N_fft,N_paths,N_bs);
Us2=zeros(M,N_paths,N_bs);
J_tau=zeros(length(Rtau_e),1);
tic
for nbs=1:N_bs
    RY1=zeros(N_fft);
    for q=1:Q
        Y0=Y_t0(:,:,q,nbs).';
        RY1=RY1+1/(Q*M)*(Y0*Y0');
    end
    [U1,tezhengzhi1,~]=svd(RY1);
    
    Us1(:,:,nbs)=U1(:,1:N_paths);
    
    for kk1=1:length(Rtau_e)
        J_tau(kk1)=tau_music(Us1(:,:,nbs),Rtau_e(kk1),N_fft,delta_f);
    end
end
toc



disp('Smusic-theta')
Rtheta_e=0.5:0.5:90;
J_theta=zeros(length(Rtheta_e),1);
tic
for nbs=1:N_bs
    RY1=zeros(N_fft);
    RY2=zeros(M);
    for q=1:Q
        RY2=RY2+1/(Q*N_fft)*(Y_t0(:,:,q,nbs)*Y_t0(:,:,q,nbs)');
    end
    
    [U2,tezhengzhi2,~]=svd(RY2);
    
    Us2(:,:,nbs)=U2(:,1:N_paths);
    
    for kk2=1:length(Rtheta_e)
        J_theta(kk2)=theta_music(Us2(:,:,nbs),B,Rtheta_e(kk2)/180*pi,N_v,N_h);
    end
end
toc

disp('spotfi-FT')
J_spotfiFT=zeros(length(Rtheta_e),length(Rtau_e));
tic
for nbs=1:N_bs
    for kk1=1:length(Rtau_e)
        for kk2=1:length(Rtheta_e)
            J_spotfiFT(kk2,kk1)=spotfi_FT0(Ydpd_t(:,:,nbs),B,[Rtau_e(kk1)/1e9;Rtheta_e(kk2)/180*pi],N_v,N_h,N_fft,delta_f);
        end
    end
end
toc

disp('spotfi-music')
J_spotfiMusic=zeros(length(Rtheta_e),length(Rtau_e));
Us=zeros(LenFreWindow*M,N_paths,N_bs);%通道输出频域数据
tic
for nbs=1:N_bs
    Y_ph=zeros(LenFreWindow*M,Qf*Q);
    Y0=Ydpd_t(:,:,nbs);
    for qf=1:Qf
        Y_ph(:,(qf-1)*Q+(1:Q))=Y0((qf-1)*M+(1:LenFreWindow*M),:);
    end
    [U,~,~]=svd(Y_ph);
    Us(:,:,nbs)=U(:,1:N_paths);
    for kk1=1:length(Rtau_e)
        for kk2=1:length(Rtheta_e)
            J_spotfiMusic(kk2,kk1)=spotfi_bs0(Us(:,:,nbs),B,[Rtau_e(kk1)/1e9;Rtheta_e(kk2)/180*pi],N_v,N_h,LenFreWindow,delta_f);
        end
    end
end
toc

xx=.1:.1:10;
yy=.1:.1:10;

J_MLdpd=zeros(length(yy),length(xx));
disp('dpd-FT')
tic
for kk1=1:length(xx)
    for kk2=1:length(yy)
        J_MLdpd(kk2,kk1)=-toaMLdpd(Ydpd_t,X0,B,P_bs,[xx(kk1),yy(kk2)],D,fn0,lamda);
    end
end
toc

J_Musicdpd=zeros(length(yy),length(xx));
disp('dpd-music')

tic
for nbs=1:N_bs
    Y_ph=zeros(LenFreWindow*M,Qf*Q);
    Y0=Ydpd_t(:,:,nbs);
    for qf=1:Qf
        Y_ph(:,(qf-1)*Q+(1:Q))=Y0((qf-1)*M+(1:LenFreWindow*M),:);
    end
    [U,~,~]=svd(Y_ph);
    Us(:,:,nbs)=U(:,1:N_paths);
end
for kk1=1:length(xx)
    for kk2=1:length(yy)
        J_Musicdpd(kk2,kk1)=-Music_dpd(Us,X0(1:LenFreWindow),B,P_bs,[xx(kk1),yy(kk2)],D,fn0(1:LenFreWindow),lamda);
    end
end
toc