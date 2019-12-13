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
N_mp=6;%�ྶ���������׾���
Res_aoa=1;%�ྶ�Ƕȼ��
Res_tau=4e-9;%ʱ�Ӽ��

delta_f=120e3*8;%OFDMƵ�ʼ��
N_fft=2048/8;%FFT����
Bw=200e6;%�źŴ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_idxRange=90/Res_aoa;
Res_aoa=Res_aoa/180*pi;

tau_idxRange=ceil((N_mp-1)/2*3);%���ʱ�ӵ�idx

SNR=20;%�����dB��%ĿǰSNR������symbol�����������������֮��
sigma2=10^(-174/10-3)*Bw;%1;%�������ʣ�W
Q=1;%������

%����������
CFO=0;%�ز�Ƶ��ƫ��1.5e3
STO=0;%������ʱƫ��
SFO=0;%����Ƶ��Ƶ��µĲ�������ƫ��0.6e-9/10e-3

T=1/delta_f;%�۲�ʱ��s,����Ƶ�ʷֱ���
fs=N_fft*delta_f;%������������ʱ��ֱ���
Ts=1/fs;%�źŲ�������
N_ofdm=N_fft;%OFDM���ز���
fn=(1:N_ofdm)*delta_f;%�����ź�Ƶ��
tt=(0:N_fft-1)*Ts;%BSʱ���ᣬ��λ��
X=sqrt(10^(SNR/10)*sigma2)*ones(N_ofdm,Q);%�ŵ���������   %exp(1j*2*pi*rand(N_ofdm,Q))�����Ƶ����;
B=DFTcodebook_spotfimuci_2(N_v,N_h,M);%eye(M);��ͨ�����������뱾
%% ���ɽ����ź�
xx=0.5:1.5:29.5;%
yy=0.5:0.5:9.5;%��������
N_simu=30;
error0=zeros(length(yy),length(xx),N_simu);
RMSE=zeros(length(yy),length(xx));
for kk1=1:length(xx)
    for kk2=1:length(yy)
        %%%%%%%%%%%%%%%%%%%%%%%%%�ź�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_ue=[xx(kk1),yy(kk2),4];%UEλ��
        gama=atan2(P_ue(2)-P_bs(:,2),P_ue(1)-P_bs(:,1));%Ŀ�굽BS�ľ�������Ƕ�
        parfor simu_idx=1:N_simu
            P_bs0=P_bs;
            gama0=gama;
            gama_H=gama_H0;
            D=D1;
            X0=X;
            fn0=fn;
            
            Ydpd_t=zeros(N_ofdm*M,Q,N_bs);%ͨ�����Ƶ������
            tau_idx=zeros(N_mp,N_bs);
            tau0=zeros(N_mp,N_bs);%�ྶ�ź�ʱ��
            theta_mp=zeros(N_mp-1,N_bs);
            theta=zeros(N_mp,N_bs);%�ྶ��Է�λ��
            %         phi=zeros(N_mp-1,N_bs);%�ྶ�߶Ƚ�
            F_alpha=ones(N_mp,N_bs);%�ྶϵ������˥������,��ʼ��Ϊ��˥��
            A=zeros(M,N_mp,N_bs);%�ྶ��������
            for nbs=1:N_bs
                tau_idx(:,nbs)=[0 1 randperm(tau_idxRange,N_mp-2)+1];
                tau_idx(2:end,nbs)=tau_idx(randperm(N_mp-1)+1,nbs);
                tau0(:,nbs)=tau_idx(:,nbs)*Res_tau+norm(P_ue-P_bs0(nbs,:))/c;
                
                theta_mp(:,nbs)=[gama0(nbs)-gama_H(nbs)+Res_aoa (randperm(theta_idxRange,N_mp-2)-(theta_idxRange/2))*Res_aoa];
                theta_mp(:,nbs)=theta_mp(randperm(N_mp-1),nbs);
                theta(:,nbs)=[gama0(nbs)-gama_H(nbs);theta_mp(:,nbs)];
                
                %             Kvec=-2*pi/lamda*[(P_ue-P_bs(nbs,:))/norm(P_ue-P_bs(nbs,:));(cos(phi(:,nbs)).*cos(gama_H(nbs)+theta(2:end,nbs))) (cos(phi(:,nbs)).*sin(gama_H(nbs)+theta(2:end,nbs))) sin(phi(:,nbs))];
                Kvec=-2*pi/lamda*[(P_ue-P_bs0(nbs,:))/norm(P_ue-P_bs0(nbs,:));(cos(gama_H(nbs)+theta(2:end,nbs))) (sin(gama_H(nbs)+theta(2:end,nbs))) zeros(N_mp-1,1)];
                A(:,:,nbs)=exp(1j*(Kvec*D(:,:,nbs)).');
                F_alpha(:,nbs)=[1,0.5+rand(1,N_mp-1)];
                alpha=diag(F_alpha(:,nbs))*exp(1j*rand(N_mp,1)*2*pi);
                
                %�����ض������½����ź��е�ȷ������
                
                for k=1:N_ofdm
                    m_t0=zeros(M,1);
                    for nmp=1:N_mp
                        m_t0=m_t0+alpha(nmp)*A(:,nmp,nbs)*exp(-1j*2*pi*fn0(k)*tau0(nmp,nbs));                       
                    end
                    Ydpd_t((k-1)*M+1:k*M,:,nbs)=B*m_t0*X0(k,:);
                end
                Ydpd_t(:,:,nbs)=Ydpd_t(:,:,nbs)+wgn(N_ofdm*M,Q,sigma2,'linear','complex');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%λ�ù���%%%%%%%%%%%%%%%%%%%%%%%%%%
%             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
            [P_ue_e,fval,~]=fminsearch(@(P_ue2)FastHR_dpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),P_ue(1:2));
            error0(simu_idx,kk2,kk1)=norm([P_ue_e,4]-P_ue);
        end
        RMSE(kk2,kk1)=sqrt(1/N_simu*sum(error0(:,kk2,kk1).^2));
    end
end
toc
  save(['..\mat\' 'FastHRdpd_RMSE_F0.5-1.5' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'Q=' num2str(Q) 'M' num2str(M*N_fft) 'Nfft.mat'],'RMSE');

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
