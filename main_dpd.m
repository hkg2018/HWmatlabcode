% clear all;
tic
%% 仿真参数
c=3e8;%光速
lamda=c/30e9;%载波波长

P_bs=[
    0,0,4;
    30,0,4;
    15,10,4
    ];%BS位置
N_bs=size(P_bs,1);%BS数量

N_v=8;%面阵行数
N_h=8;%面阵列数
d=lamda/2;%阵元间距
M=N_h*N_v;%阵子数
gama_H0=[pi/4 3*pi/4 -pi/2];%各个天线阵列的旋转角度=绝对入射角度-相对入射角度

D0=zeros(3,M);%各阵元的相对位置坐标(BS阵列形状)
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
        ]*D0;%旋转天线阵列
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_mp=6;%多径数量（含首径）
Res_aoa=1;%多径角度间隔
Res_tau=4e-9;%时延间隔

delta_f=120e3*8;%OFDM频率间隔
N_fft=2048/8;%FFT点数
Bw=200e6;%信号带宽
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_idxRange=90/Res_aoa;
Res_aoa=Res_aoa/180*pi;

tau_idxRange=ceil((N_mp-1)/2*3);%最大时延的idx

SNR=20;%信噪比dB，%目前SNR定义是symbol功率与测量噪声功率之比
sigma2=10^(-174/10-3)*Bw;%1;%噪声功率，W
Q=1;%样本数

%非理想因素
CFO=0;%载波频率偏差1.5e3
STO=0;%采样定时偏差
SFO=0;%采样频率频差导致的采样周期偏差0.6e-9/10e-3

T=1/delta_f;%观测时长s,决定频率分辨率
fs=N_fft*delta_f;%基带带宽，决定时间分辨率
Ts=1/fs;%信号采样周期
N_ofdm=N_fft;%OFDM子载波数
fn=(1:N_ofdm)*delta_f;%基带信号频点
tt=(0:N_fft-1)*Ts;%BS时间轴，单位秒
X=sqrt(10^(SNR/10)*sigma2)*ones(N_ofdm,Q);%信道测试序列   %exp(1j*2*pi*rand(N_ofdm,Q))随机导频序列;
B=DFTcodebook_spotfimuci_2(N_v,N_h,M);%eye(M);单通道正交波束码本
%% 生成接收信号
xx=0.5:1.5:29.5;%
yy=0.5:0.5:9.5;%场景网格
N_simu=30;
error0=zeros(length(yy),length(xx),N_simu);
RMSE=zeros(length(yy),length(xx));
for kk1=1:length(xx)
    for kk2=1:length(yy)
        %%%%%%%%%%%%%%%%%%%%%%%%%信号生成%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_ue=[xx(kk1),yy(kk2),4];%UE位置
        gama=atan2(P_ue(2)-P_bs(:,2),P_ue(1)-P_bs(:,1));%目标到BS的绝对入射角度
        parfor simu_idx=1:N_simu
            P_bs0=P_bs;
            gama0=gama;
            gama_H=gama_H0;
            D=D1;
            X0=X;
            fn0=fn;
            
            Ydpd_t=zeros(N_ofdm*M,Q,N_bs);%通道输出频域数据
            tau_idx=zeros(N_mp,N_bs);
            tau0=zeros(N_mp,N_bs);%多径信号时延
            theta_mp=zeros(N_mp-1,N_bs);
            theta=zeros(N_mp,N_bs);%多径相对方位角
            %         phi=zeros(N_mp-1,N_bs);%多径高度角
            F_alpha=ones(N_mp,N_bs);%多径系数方差衰减因子,初始化为不衰减
            A=zeros(M,N_mp,N_bs);%多径阵列流形
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
                
                %生成特定场景下接收信号中的确定部分
                
                for k=1:N_ofdm
                    m_t0=zeros(M,1);
                    for nmp=1:N_mp
                        m_t0=m_t0+alpha(nmp)*A(:,nmp,nbs)*exp(-1j*2*pi*fn0(k)*tau0(nmp,nbs));                       
                    end
                    Ydpd_t((k-1)*M+1:k*M,:,nbs)=B*m_t0*X0(k,:);
                end
                Ydpd_t(:,:,nbs)=Ydpd_t(:,:,nbs)+wgn(N_ofdm*M,Q,sigma2,'linear','complex');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位置估计%%%%%%%%%%%%%%%%%%%%%%%%%%
%             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
            [P_ue_e,fval,~]=fminsearch(@(P_ue2)FastHR_dpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),P_ue(1:2));
            error0(simu_idx,kk2,kk1)=norm([P_ue_e,4]-P_ue);
        end
        RMSE(kk2,kk1)=sqrt(1/N_simu*sum(error0(:,kk2,kk1).^2));
    end
end
toc
  save(['..\mat\' 'FastHRdpd_RMSE_F0.5-1.5' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'Q=' num2str(Q) 'M' num2str(M*N_fft) 'Nfft.mat'],'RMSE');

%% 绘制位置谱
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
