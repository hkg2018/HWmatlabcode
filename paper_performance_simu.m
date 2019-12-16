% clear all;
tic
%% �������
c=3e8;%����
lamda=c/30e9;%�ز�����

P_bs=[
    0,0,4;
    30,0,4;
    15,10,4
    ];%BSλ��
N_bs=size(P_bs,1);%BS����

N_v=8;%��������
N_h=8;%��������
d=lamda/2;%��Ԫ���
M=N_h*N_v;%������
gama_H0=[pi/4 3*pi/4 -pi/2];%�����������е���ת�Ƕ�=��������Ƕ�-�������Ƕ�

D0=zeros(3,M);%����Ԫ�����λ������(BS������״)
D1=zeros(3,M,N_bs);
for nh=1:N_h
    for nv=1:N_v
        D0(:,(nh-1)*N_v+nv)=[(nv-1)*d;(nh-1)*d;0];
    end
end
for nbs=1:N_bs
    D1(:,:,nbs)=[
        cos(gama_H0(nbs)),sin(-gama_H0(nbs)),0;
        -sin(-gama_H0(nbs)),cos(gama_H0(nbs)),0;
        0,0,1
        ]*D0;%��ת��������
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_paths=6;%�ྶ���������׾���
Res_aoa=2;%�ྶ�Ƕȼ��
Res_tau0=[2,5]*1e-9;%ʱ�Ӽ��

delta_f=120e3*8;%OFDMƵ�ʼ��
N_fft=2048/8;%FFT����
Bw=200e6;%�źŴ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_idxRange=90/Res_aoa;
Res_aoa=Res_aoa/180*pi;

tau_idxRange=ceil((N_paths-1)/2*3);%���ʱ�ӵ�idx

SNR=[0,20];%�����dB��%ĿǰSNR������symbol�����������������֮��
sigma2=10^(-174/10-3)*Bw;%1;%�������ʣ�W
Q=1;%������

%����������
CFO0=0;%��һ���ز�Ƶ��ƫ��
% CFO_org=CFO0(ii);
CFO=zeros(1,M);
STO=0;%������ʱƫ��
SFO=0;%����Ƶ��Ƶ��µĲ�������ƫ��0.6e-9/10e-3

T=1/delta_f;%�۲�ʱ��s,����Ƶ�ʷֱ���
fs=N_fft*delta_f;%������������ʱ��ֱ���
Ts=1/fs;%�źŲ�������
N_ofdm=N_fft;%OFDM���ز���
fn=(1:N_ofdm)*delta_f;%�����ź�Ƶ��
tt=(0:N_fft-1)*Ts;%BSʱ���ᣬ��λ��
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);��ͨ�����������뱾
%% ���ɽ����ź�
xx=0.5:1.5:29.5;%
yy=0.5:0.5:9.5;%��������
N_simu=1;
LenFreWindow=fix(N_fft*0.25);
Qf=N_fft-LenFreWindow+1;
MFerror0=zeros(N_simu,length(yy),length(xx));
MusiCerror0=zeros(N_simu,length(yy),length(xx));
MFrmse=zeros(length(yy),length(xx));
MusiCrmse=zeros(length(yy),length(xx));
Pos_err2=zeros(5,N_simu,length(yy),length(xx));
Pos_rmse=zeros(5,length(yy),length(xx));
for idx=1:length(SNR)
    X=sqrt(10^(SNR(idx)/10)*sigma2)*ones(N_ofdm,Q);%�ŵ���������   %exp(1j*2*pi*rand(N_ofdm,Q))�����Ƶ����;
    for ii=1:length(Res_tau0)
        Res_tau=Res_tau0(ii);       
        for simu_idx=1:N_simu
            for kk1=1:length(xx)
                xx0=xx(kk1);
                parfor kk2=1:length(yy)

                    P_bs0=P_bs;
                    gama_H=gama_H0;
                    D=D1;
                    X0=X;
                    fn0=fn;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%�ź�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    P_ue=[xx0,yy(kk2),4];%UEλ��
                    gama0=atan2(P_ue(2)-P_bs0(:,2),P_ue(1)-P_bs0(:,1));%Ŀ�굽BS�ľ�������Ƕ�
                    
                    Ydpd_t=zeros(N_ofdm*M,Q,N_bs);%ͨ�����Ƶ������
                    Y_t=zeros(M,Q,N_ofdm,N_bs);%ͨ�����Ƶ������
                    Us=zeros(LenFreWindow*M,N_paths,N_bs);%ͨ�����Ƶ������
                    tau_idx=zeros(N_paths,N_bs);
                    tau0=zeros(N_paths,N_bs);%�ྶ�ź�ʱ��
                    theta_mp=zeros(N_paths-1,N_bs);
                    theta=zeros(N_paths,N_bs);%�ྶ��Է�λ��
                    %         phi=zeros(N_mp-1,N_bs);%�ྶ�߶Ƚ�
                    F_alpha=ones(N_paths,N_bs);%�ྶϵ������˥������,��ʼ��Ϊ��˥��
                    A=zeros(M,N_paths,N_bs);%�ྶ��������
                    for nbs=1:N_bs
                        tau_idx(:,nbs)=[0 1 randperm(tau_idxRange,N_paths-2)+1];
%                         tau_idx(2:end,nbs)=tau_idx(randperm(N_paths-1)+1,nbs);
                        tau0(:,nbs)=diag([0 Res_tau ones(1,N_paths-2)*1e-9])*tau_idx(:,nbs)+norm(P_ue-P_bs0(nbs,:))/c;
                        
                        theta_mp(:,nbs)=[gama0(nbs)-gama_H(nbs)+Res_aoa (randperm(theta_idxRange,N_paths-2)-(theta_idxRange/2))*Res_aoa];
                        theta_mp(:,nbs)=theta_mp(randperm(N_paths-1),nbs);
                        theta(:,nbs)=[gama0(nbs)-gama_H(nbs);theta_mp(:,nbs)];
                        
                        %             Kvec=-2*pi/lamda*[(P_ue-P_bs(nbs,:))/norm(P_ue-P_bs(nbs,:));(cos(phi(:,nbs)).*cos(gama_H(nbs)+theta(2:end,nbs))) (cos(phi(:,nbs)).*sin(gama_H(nbs)+theta(2:end,nbs))) sin(phi(:,nbs))];
                        Kvec=-2*pi/lamda*[(P_ue-P_bs0(nbs,:))/norm(P_ue-P_bs0(nbs,:));(cos(gama_H(nbs)+theta(2:end,nbs))) (sin(gama_H(nbs)+theta(2:end,nbs))) zeros(N_paths-1,1)];
                        A(:,:,nbs)=exp(1j*(Kvec*D(:,:,nbs)).');
                        F_alpha(:,nbs)=[1,0.5+0.6*rand(1,N_paths-1)];
                        alpha=diag(F_alpha(:,nbs))*exp(1j*rand(N_paths,1)*2*pi);
%                         if CFO_org
%                             CFO=CFO_org*(2*rand(1,M)-1);
%                         else
%                             CFO=0;
%                         end
                        
                        %�����ض������½����ź��е�ȷ������
                        
                        for k=1:N_ofdm
                            H0=zeros(M,1);
                            for nmp=1:N_paths
                                H0=H0+alpha(nmp)*A(:,nmp,nbs)*exp(-1j*2*pi*fn0(k)*tau0(nmp,nbs));
                            end
                            Yt0=B*H0*X0(k,:);%diag(exp(1j*2*pi*fix((k-1)/8)*CFO))*
                            Ydpd_t((k-1)*M+1:k*M,:,nbs)=Yt0;
                            Y_t(:,:,k,nbs)=Yt0+wgn(M,Q,sigma2,'linear','complex');
                        end
                        Ydpd_t(:,:,nbs)=Ydpd_t(:,:,nbs)+wgn(N_ofdm*M,Q,sigma2,'linear','complex');
%                         if Qf~=1
                            Y_ph=zeros(LenFreWindow*M,Qf*Q);
                            Y0=Ydpd_t(:,:,nbs);
                            for qf=1:Qf
                                Y_ph(:,(qf-1)*Q+(1:Q))=Y0((qf-1)*M+(1:LenFreWindow*M),:);
                            end
                            [U,~,~]=svd(Y_ph);
                            Us(:,:,nbs)=U(:,1:N_paths);
%                         end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%λ�ù���%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
%                     if Qf~=1
                        %Musicdpd
                        [MusicP_ue,~,~]=fminsearch(@(P_ue2)Music_dpd(Us,X0(1:LenFreWindow),B,P_bs0,P_ue2,D,fn0(1:LenFreWindow),lamda),P_ue(1:2));
                        MusiCerror0(simu_idx,kk2,kk1)=norm([MusicP_ue,4]-P_ue);
                        
                        %spotfi_Music
                        [spotfi_tau,spotfi_Theta]=spotfi_estimator(Us,B,N_v,N_h,delta_f,tau0,theta,N_paths,gama_H);
%                     else
                        %FastHRdpd
                        [FastHRP_ue,~,~]=fminsearch(@(P_ue2)FastHR_dpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),P_ue(1:2));
                        [MFP_ue,~,~]=fminsearch(@(P_ue2)MFdpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),FastHRP_ue);
                        FastHRerror0(simu_idx,kk2,kk1)=norm([FastHRP_ue,4]-P_ue);
                        MFerror0(simu_idx,kk2,kk1)=norm([MFP_ue,4]-P_ue);
                        %spotfi_MF
                        [spotfiMF_tau,spotfiMF_Theta]=spotfi_MFestimator(Ydpd_t,B,N_v,N_h,N_fft,delta_f,tau0,theta,N_paths,gama_H);
%                     end
                    
                    %tau_music & theta_music
                    [MUSIC_tau,MUSIC_Theta]=tauTheta_MUSICestimator(Y_t,B,N_v,N_h,N_fft,delta_f,tau0,theta,N_paths,gama_H);
                    %����abs(([spotfi_tau;MUSIC_tau]-tau0(1,:))*1e9);abs([spotfi_Theta;MUSIC_Theta]-theta(1,:)/pi*180)
                    Pos_err2(:,simu_idx,kk2,kk1)=spotfi_Location_joint1203_1([spotfi_Theta;spotfi_tau;MUSIC_Theta;MUSIC_tau;spotfiMF_Theta;spotfiMF_tau],P_bs0,P_ue);               
                end
            end
        end
        MusiCerror=permute(MusiCerror0,[2,3,1]);
        MusiCrmse=sqrt(1/N_simu*sum(MusiCerror.^2,3));
        FastHRerror=permute(FastHRerror0,[2,3,1]);
        FastHRrmse=sqrt(1/N_simu*sum(FastHRerror.^2,3));
        MFerror=permute(MFerror0,[2,3,1]);
        MFrmse=sqrt(1/N_simu*sum(MFerror.^2,3));
        Pos_2err=permute(Pos_err2,[3,4,1,2]);
        Pos_rmse=sqrt(sum(Pos_2err,4)/N_simu);
        save(['..\mat\' 'performance_all_estimator' num2str(Res_aoa/pi*180) 'deg'  num2str(Res_tau*1e9) 'ns_' 'SNR' num2str(SNR(idx)) 'dB.mat'],'MusiCrmse','FastHRrmse','MFrmse','Pos_rmse');
    end
end
toc
for idx=1:length(Res_tau0)
    figure(idx);
    y=zeros(length(SNR),8);
    for ii=1:length(SNR)
        load(['..\mat\' 'performance_all_estimator' num2str(Res_aoa/pi*180) 'deg'  num2str(Res_tau0(idx)*1e9) 'ns_' 'SNR' num2str(SNR(ii)) 'dB.mat']);
        
        OtherRmse=zeros(size(Pos_rmse,1)*size(Pos_rmse,2),5);
        for jj=1:size(OtherRmse,2)
            rmse0=Pos_rmse(:,:,jj);
            OtherRmse(:,jj)=rmse0(:);
        end
        [x,f]=ploTcdf([MusiCrmse(:),FastHRrmse(:),OtherRmse,MFrmse(:)]);
        
        x=x(:,[1,3,2,8,6,7,4,5]);
        f=f(:,[1,3,2,8,6,7,4,5]);
        for jj=1:size(y,2)
        y(ii,jj)=x(find(f(:,jj)>=0.9,1),jj);
        end       
    end
        semilogy(SNR.'*ones(1,size(y,2)),y,'LineWidth',2);
        grid on;
        xlabel('SNR/dB');
        ylabel('error(m)@90%');
%         title({[ 'normalizedCFO=' num2str(CFO0(ii)) '(' num2str(CFO0(ii)*delta_f) 'Hz)'];[ 'ResTau=' num2str(Res_tau0(ii)*1e9) 'ns,' 'ResAoA=' num2str(Res_aoa/pi*180) 'deg']});
        legend('dpd-Music','Spotfi-Music','iprVdpd-MF','dpd-MF','TOA-solution','Spotfi-MF','Joint-solution','AOA-solution(MaxPower)','location','best');
%         axis([x(1,1),0.3,0,100]);
%         set(gca,'xtick',[.1,.2,.3]);
end


%% ����λ����
% tic
% xx=1:1:29;
% yy=1:1:29;
%     J_FastHRdpd=zeros(length(yy),length(xx));
%     J_MLdpd=zeros(length(yy),length(xx));
%     for kk1=1:length(xx)
%         for kk2=1:length(yy)
%             J_MLdpd(kk2,kk1)=-toaMLdpd(Ydpd_t,X0,B,P_bs,[xx(kk1),yy(kk2)],D,fn0,lamda);
%             J_FastHRdpd(kk2,kk1)=-FastHR_dpd(Ydpd_t,X0,B,P_bs,[xx(kk1),yy(kk2)],D,fn0,lamda);
%         end
%     end
%     J_MLdpd=J_MLdpd/max(max(J_MLdpd));
%     J_FastHRdpd=J_FastHRdpd/max(max(J_FastHRdpd));
%
%
%     [XmirrorUE,YmirrorUE]=pol2cart(theta+gama_H,tau0*c);
%     XmirrorUE=XmirrorUE+P_bs(:,1).';
%     YmirrorUE=YmirrorUE+P_bs(:,2).';
%
%     fig(1)=figure;
%     contour(xx,yy,J_MLdpd);
%     hold on;
%     scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
%     scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
%     scatter3(P_ue_e(1),P_ue_e(2),1,100,'b','+');
%     scatter3(XmirrorUE(2:end,1),YmirrorUE(2:end,1),ones(N_mp-1,1),50,'r','p');
%     scatter3(XmirrorUE(2:end,2),YmirrorUE(2:end,2),ones(N_mp-1,1),50,'b','p');
%     scatter3(XmirrorUE(2:end,3),YmirrorUE(2:end,3),ones(N_mp-1,1),50,'g','p');
%     grid on;
%     title(['MLdpd position spectrum' 'N_{paths}=' num2str(N_mp) ',Tau-Theta_{res}[' num2str(Res_tau*1e9) ',' num2str(Res_aoa/pi*180) '],'...
%           'error=' num2str(error) 'm']);
%     xlabel('x/m');
%     ylabel('y/m');
%     axis([0,30,0,30]);
%     legend('Power','BS','UE','estimated position', '1#mirrorUE','2#mirrorUE','3#mirrorUE', 'Location','best');
%
%     fig(2)=figure;
%     contour(xx,yy,J_FastHRdpd);
%     hold on;
%     scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
%     scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
%     scatter3(P_ue_e(1),P_ue_e(2),1,100,'b','+');
%     scatter3(XmirrorUE(2:end,1),YmirrorUE(2:end,1),ones(N_mp-1,1),50,'r','p');
%     scatter3(XmirrorUE(2:end,2),YmirrorUE(2:end,2),ones(N_mp-1,1),50,'b','p');
%     scatter3(XmirrorUE(2:end,3),YmirrorUE(2:end,3),ones(N_mp-1,1),50,'g','p');
%     grid on;
%     title({'FastHRdpd position spectrum';[ 'N_{paths}=' num2str(N_mp) ',Tau-Theta_{res}[' num2str(Res_tau*1e9) ',' num2str(Res_aoa/pi*180) '],'...
%           'error=' num2str(error) 'm']});
%     xlabel('x/m');
%     ylabel('y/m');
%     axis([0,30,0,30]);
%     legend('Power','BS','UE','estimated position', '1#mirrorUE','2#mirrorUE','3#mirrorUE', 'Location','best');
%
% toc
