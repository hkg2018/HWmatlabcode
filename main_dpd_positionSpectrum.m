% clear all;
tic
%% 仿真参数
c=3e8;%光速
lamda=c/30e9;%载波波长

P_bs=[
    0,0,4;
    30,0,4;
    15,30,4
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
N_paths=6;%多径数量（含首径）
Res_aoa=1;%多径角度间隔
Res_tau=1e-9;%时延间隔

delta_f=120e3*8;%OFDM频率间隔
N_fft=2048/8;%FFT点数
Bw=200e6;%信号带宽
Q=1;%样本数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_idxRange=90/Res_aoa;
Res_aoa=Res_aoa/180*pi;

tau_idxRange=ceil((N_paths-1)/2*3);%最大时延的idx

SNR=20;%信噪比dB，%目前SNR定义是symbol功率与测量噪声功率之比
sigma2=10^(-174/10-3)*Bw;%1;%噪声功率，W

%非理想因素
CFO0=0;%载波频率偏差Hz
% N_channel=1;
CFO=zeros(1,M);
STO=0;%采样定时偏差
SFO=0;%采样频率频差导致的采样周期偏差0.6e-9/10e-3

T=1/delta_f;%观测时长s,决定频率分辨率
fs=N_fft*delta_f;%基带带宽，决定时间分辨率
Ts=1/fs;%信号采样周期
N_ofdm=N_fft;%OFDM子载波数
fn=(1:N_ofdm)*delta_f;%基带信号频点
tt=(0:N_fft-1)*Ts;%BS时间轴，单位秒
X=sqrt(10^(SNR/10)*sigma2)*ones(N_ofdm,Q);%信道测试序列   %exp(1j*2*pi*rand(N_ofdm,Q))随机导频序列;
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);单通道正交波束码本
%% 生成接收信号
xx=15;%0.5:1.5:29.5
yy=15;%0.5:0.5:9.5场景网格
N_simu=1;
LenFreWindow=fix(N_fft*0.25);
Qf=N_fft-LenFreWindow+1;
error1=zeros(length(yy),length(xx),N_simu);
RMSE1=zeros(length(yy),length(xx));
error2=zeros(length(yy),length(xx),N_simu);
RMSE2=zeros(length(yy),length(xx));
error3=zeros(length(yy),length(xx),N_simu);
RMSE3=zeros(length(yy),length(xx));
for kk1=1:length(xx)
    for kk2=1:length(yy)
        %信号生成%
        P_ue=[xx(kk1),yy(kk2),4];%UE位置
        gama=atan2(P_ue(2)-P_bs(:,2),P_ue(1)-P_bs(:,1));%目标到BS的绝对入射角度
        for simu_idx=1:N_simu
            P_bs0=P_bs;
            gama0=gama;
            gama_H=gama_H0;
            D=D1;
            X0=X;
            fn0=fn;
            
            Ydpd_t=zeros(N_ofdm*M,Q,N_bs);%通道输出频域数据
            Y_t=zeros(M,Q,N_ofdm,N_bs);%通道输出频域数据
            if Qf~=1
            Us=zeros(LenFreWindow*M,N_paths,N_bs);%通道输出频域数据
            end
            tau_idx=zeros(N_paths,N_bs);
            tau0=zeros(N_paths,N_bs);%多径信号时延
            theta_mp=zeros(N_paths-1,N_bs);
            theta=zeros(N_paths,N_bs);%多径相对方位角
            %         phi=zeros(N_mp-1,N_bs);%多径高度角
            F_alpha=ones(N_paths,N_bs);%多径系数方差衰减因子,初始化为不衰减
            A=zeros(M,N_paths,N_bs);%多径阵列流形

            for nbs=1:N_bs
                %                 tau_idx(2:end,nbs)=tau_idx(randperm(N_mp-1)+1,nbs);
                tau0(:,nbs)=[0 Res_tau 1e-9*randperm(tau_idxRange,N_paths-2)+Res_tau]+norm(P_ue-P_bs0(nbs,:))/c;
                
                theta_mp(:,nbs)=[gama0(nbs)-gama_H(nbs)+Res_aoa (randperm(theta_idxRange,N_paths-2)-(theta_idxRange/2))*Res_aoa];
                theta_mp(:,nbs)=theta_mp(randperm(N_paths-1),nbs);
                theta(:,nbs)=[gama0(nbs)-gama_H(nbs);theta_mp(:,nbs)];
%                 theta/pi*180
                
                %             Kvec=-2*pi/lamda*[(P_ue-P_bs(nbs,:))/norm(P_ue-P_bs(nbs,:));(cos(phi(:,nbs)).*cos(gama_H(nbs)+theta(2:end,nbs))) (cos(phi(:,nbs)).*sin(gama_H(nbs)+theta(2:end,nbs))) sin(phi(:,nbs))];
                Kvec=-2*pi/lamda*[(P_ue-P_bs0(nbs,:))/norm(P_ue-P_bs0(nbs,:));(cos(gama_H(nbs)+theta(2:end,nbs))) (sin(gama_H(nbs)+theta(2:end,nbs))) zeros(N_paths-1,1)];
                A(:,:,nbs)=exp(1j*(Kvec*D(:,:,nbs)).');
                F_alpha(:,nbs)=[1,0.6+0.6*rand(1,N_paths-1)];
                alpha=diag(F_alpha(:,nbs))*exp(1j*rand(N_paths,1)*2*pi);
                if CFO0
                    CFO=CFO0*(2*rand(1,M)-1); 
                else
                    CFO=0;
                end
                
                %生成特定场景下接收信号中的确定部分
                
                for k=1:N_ofdm
                    H0=zeros(M,1);%信道空-频域响应
                    for nmp=1:N_paths
                        H0=H0+alpha(nmp)*A(:,nmp,nbs)*exp(-1j*2*pi*fn0(k)*tau0(nmp,nbs));
                    end
                    Yt0=diag(exp(1j*2*pi*fix((k-1)/8)*CFO/delta_f))*B*H0*X0(k,:);
                    Ydpd_t((k-1)*M+1:k*M,:,nbs)=Yt0;
                    Y_t(:,:,k,nbs)=Yt0+wgn(M,Q,sigma2,'linear','complex');
                end
                Ydpd_t(:,:,nbs)=Ydpd_t(:,:,nbs)+wgn(N_ofdm*M,Q,sigma2,'linear','complex');
                if Qf~=1
                    Y_ph=zeros(LenFreWindow*M,Qf*Q);
                    Y0=Ydpd_t(:,:,nbs);
                    for qf=1:Qf
                        Y_ph(:,(qf-1)*Q+(1:Q))=Y0((qf-1)*M+(1:LenFreWindow*M),:);
                    end
                    [U,~,~]=svd(Y_ph);
                    Us(:,:,nbs)=U(:,1:N_paths);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位置估计%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
            if Qf~=1
                [MusicP_ue,~,~]=fminsearch(@(P_ue2)Music_dpd(Us,X0(1:LenFreWindow),B,P_bs0,P_ue2,D,fn0(1:LenFreWindow),lamda),P_ue(1:2));
                error3(simu_idx,kk2,kk1)=norm([MusicP_ue,4]-P_ue);
            end
            [FastHRP_ue,~,~]=fminsearch(@(P_ue2)FastHR_dpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),P_ue(1:2));
            [MLP_ue,~,~]=fminsearch(@(P_ue2)MFdpd(Ydpd_t,X0,B,P_bs0,P_ue2,D,fn0,lamda),FastHRP_ue);
            
            error1(simu_idx,kk2,kk1)=norm([MLP_ue,4]-P_ue);
            error2(simu_idx,kk2,kk1)=norm([FastHRP_ue,4]-P_ue);
        end
        RMSE1(kk2,kk1)=sqrt(1/N_simu*sum(error1(:,kk2,kk1).^2));
        RMSE2(kk2,kk1)=sqrt(1/N_simu*sum(error2(:,kk2,kk1).^2));
        RMSE3(kk2,kk1)=sqrt(1/N_simu*sum(error3(:,kk2,kk1).^2));
    end
end
toc
% save(['..\mat\' 'FastHRdpd_RMSE_F1vs1.5' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat'],'RMSE');

%% 绘制位置谱
tic
xx=1:1:29;
yy=1:1:29;
J_FastHRdpd=zeros(length(yy),length(xx));
J_MLdpd=zeros(length(yy),length(xx));

for kk1=1:length(xx)
    for kk2=1:length(yy)
        J_MLdpd(kk2,kk1)=-MFdpd(Ydpd_t,X0,B,P_bs,[xx(kk1),yy(kk2)],D,fn0,lamda);
        J_FastHRdpd(kk2,kk1)=-FastHR_dpd(Ydpd_t,X0,B,P_bs,[xx(kk1),yy(kk2)],D,fn0,lamda);
    end
end
J_MLdpd=J_MLdpd/max(max(J_MLdpd));
J_FastHRdpd=J_FastHRdpd/max(max(J_FastHRdpd));


[XmirrorUE,YmirrorUE]=pol2cart(theta+gama_H,tau0*c);
XmirrorUE=XmirrorUE+P_bs(:,1).';
YmirrorUE=YmirrorUE+P_bs(:,2).';

figure('NumberTitle', 'off', 'Name', ...
    [ 'N_{paths}=' num2str(N_paths) ',Tau-Theta_{res}[' num2str(Res_tau*1e9) 'ns,' num2str(Res_aoa/pi*180) 'deg]' 'CFO' num2str(CFO0) 'Hz']);
subplot(1,3,1);
contour(xx,yy,J_MLdpd);
hold on;
scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
scatter3(MLP_ue(1),MLP_ue(2),1,100,'b','+');

scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
scatter3(XmirrorUE(2:end,1),YmirrorUE(2:end,1),ones(N_paths-1,1),50,'r','p');
scatter3(XmirrorUE(2:end,2),YmirrorUE(2:end,2),ones(N_paths-1,1),50,'b','p');
scatter3(XmirrorUE(2:end,3),YmirrorUE(2:end,3),ones(N_paths-1,1),50,'g','p');

grid on;
title({'dpd-FT position spectrum';...
    ['error=' num2str(RMSE1) 'm']});
xlabel('x/m');
ylabel('y/m');

axis([0,30,0,30]);

subplot(1,3,2);
contour(xx,yy,J_FastHRdpd);
hold on;
scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
scatter3(FastHRP_ue(1),FastHRP_ue(2),1,100,'b','+');
scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
scatter3(XmirrorUE(2:end,1),YmirrorUE(2:end,1),ones(N_paths-1,1),50,'r','p');
scatter3(XmirrorUE(2:end,2),YmirrorUE(2:end,2),ones(N_paths-1,1),50,'b','p');
scatter3(XmirrorUE(2:end,3),YmirrorUE(2:end,3),ones(N_paths-1,1),50,'g','p');
grid on;
title({'improved dpd-FT position spectrum';...
    ['error=' num2str(RMSE2) 'm']});
xlabel('x/m');
ylabel('y/m');
axis([0,30,0,30]);



if  Qf~=1
J_Musicdpd=zeros(length(yy),length(xx));
for kk1=1:length(xx)
    for kk2=1:length(yy)        
        J_Musicdpd(kk2,kk1)=-Music_dpd(Us,X0(1:LenFreWindow),B,P_bs,[xx(kk1),yy(kk2)],D,fn0(1:LenFreWindow),lamda);
    end
end
J_Musicdpd=J_Musicdpd/max(max(J_Musicdpd));

subplot(1,3,3);
contour(xx,yy,J_Musicdpd);
hold on;
scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
scatter3(MusicP_ue(1),MusicP_ue(2),1,100,'b','+');

scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
scatter3(XmirrorUE(2:end,1),YmirrorUE(2:end,1),ones(N_paths-1,1),50,'r','p');
scatter3(XmirrorUE(2:end,2),YmirrorUE(2:end,2),ones(N_paths-1,1),50,'b','p');
scatter3(XmirrorUE(2:end,3),YmirrorUE(2:end,3),ones(N_paths-1,1),50,'g','p');

grid on;
title({'dpd-music position spectrum';...
    ['error=' num2str(RMSE3) 'm']});
xlabel('x/m');
ylabel('y/m');
legend('Power','UE','estimated position','BS', '1#mirrorUE','2#mirrorUE','3#mirrorUE', 'Location','best');
axis([0,30,0,30]);
end
toc