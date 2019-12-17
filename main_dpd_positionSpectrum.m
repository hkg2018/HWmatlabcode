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
gama_H0=[pi/4,3*pi/4,-pi/2];%各个天线阵列的旋转角度=绝对入射角度-相对入射角度

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
Power_noise=-90;%dBm
sigma2=10^(Power_noise/10-3);%1;%噪声功率，W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



delta_f=120e3*8;%OFDM频率间隔
N_fft=2048/8;%FFT点数
Bw=N_fft*delta_f;%信号带宽
Q=1;%样本数
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
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);单通道正交波束码本
LenFreWindow=fix(N_fft*0.25);
Qf=N_fft-LenFreWindow+1;
%% 生成接收信号
P_ue=[15,15,4];%UE位置
gama=atan2(P_ue(2)-P_bs(:,2),P_ue(1)-P_bs(:,1)).';%目标到BS的绝对入射角度

%%%%%%%%%%%%%%%%%%%3GPP_UMI_channel_model%%%%%%%%%%%%%%%%%%%%%
theta0=gama-gama_H0;
channel_generate1217_3BS(theta0,N_bs);
% load('CSI.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ydpd_t=zeros(N_ofdm*M,Q,N_bs);%通道输出频域数据
Y_t=zeros(M,Q,N_ofdm,N_bs);%通道输出频域数据
if Qf~=1
    Us=zeros(LenFreWindow*M,max(N_paths),N_bs);%通道输出频域数据
end

for nbs=1:N_bs
    %                 theta/pi*180
    %             Kvec=-2*pi/lamda*[(P_ue-P_bs(nbs,:))/norm(P_ue-P_bs(nbs,:));(cos(phi(:,nbs)).*cos(gama_H(nbs)+theta(2:end,nbs))) (cos(phi(:,nbs)).*sin(gama_H(nbs)+theta(2:end,nbs))) sin(phi(:,nbs))];
    Kvec=-2*pi/lamda*[(P_ue-P_bs(nbs,:))/norm(P_ue-P_bs(nbs,:));(cos(gama_H0(nbs)+theta(2:N_paths(nbs),nbs))) (sin(gama_H0(nbs)+theta(2:N_paths(nbs),nbs))) zeros(N_paths(nbs)-1,1)];
    A=exp(1j*(Kvec*D1(:,:,nbs)).');%多径阵列流形
    
    if CFO0
        CFO=CFO0*(2*rand(1,M)-1);
    else
        CFO=0;
    end
    
    %生成特定场景下接收信号中的确定部分
    
    for k=1:N_ofdm
        H0=zeros(M,1);%信道空-频域响应
        for nmp=1:N_paths(nbs)
            H0=H0+alpha(nmp,nbs)*A(:,nmp)*exp(-1j*2*pi*fn(k)*tau0(nmp,nbs));
        end
        Yt0=diag(exp(1j*2*pi*fix((k-1)/8)*CFO/delta_f))*B*H0*ones(1,Q);
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
        Us(:,:,nbs)=U(:,1:N_paths(nbs));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位置估计%%%%%%%%%%%%%%%%%%%%%%%%%%
%             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
if Qf~=1
    [MusicP_ue,~,~]=fminsearch(@(P_ue2)Music_dpd(Us,N_paths,B,P_bs,P_ue2,D1,fn(1:LenFreWindow).',lamda),P_ue(1:2));
    error3=norm([MusicP_ue,4]-P_ue);
end
[FastHRP_ue,~,~]=fminsearch(@(P_ue2)FastHR_dpd(Ydpd_t,B,P_bs,P_ue2,D1,fn.',lamda),P_ue(1:2));
[MLP_ue,~,~]=fminsearch(@(P_ue2)MFdpd(Ydpd_t,B,P_bs,P_ue2,D1,fn.',lamda),FastHRP_ue);

error1=norm([MLP_ue,4]-P_ue);
error2=norm([FastHRP_ue,4]-P_ue);

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
        J_MLdpd(kk2,kk1)=-MFdpd(Ydpd_t,B,P_bs,[xx(kk1),yy(kk2)],D1,fn.',lamda);
        J_FastHRdpd(kk2,kk1)=-FastHR_dpd(Ydpd_t,B,P_bs,[xx(kk1),yy(kk2)],D1,fn.',lamda);
    end
end
J_MLdpd=J_MLdpd/max(max(J_MLdpd));
J_FastHRdpd=J_FastHRdpd/max(max(J_FastHRdpd));


[XmirrorUE,YmirrorUE]=pol2cart(theta+gama_H0,tau0*c);
XmirrorUE=XmirrorUE+P_bs(:,1).';
YmirrorUE=YmirrorUE+P_bs(:,2).';

figure('NumberTitle', 'off', 'Name', ...
    [ 'N_{paths}=' num2str(N_paths) ',Tau-Theta_{res}[' num2str(min_tau) 'ns,' num2str(min_aoa) 'deg]' 'CFO' num2str(CFO0) 'Hz']);
subplot(1,3,1);
contour(xx,yy,J_MLdpd);
hold on;
scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
scatter3(MLP_ue(1),MLP_ue(2),1,100,'b','+');

scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');

grid on;
title({'dpd-FT position spectrum';...
    ['error=' num2str(error1) 'm']});
xlabel('x/m');
ylabel('y/m');

axis([0,30,0,30]);

subplot(1,3,2);
contour(xx,yy,J_FastHRdpd);
hold on;
scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
scatter3(FastHRP_ue(1),FastHRP_ue(2),1,100,'b','+');
scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');
grid on;
title({'improved dpd-FT position spectrum';...
    ['error=' num2str(error2) 'm']});
xlabel('x/m');
ylabel('y/m');
axis([0,30,0,30]);



if  Qf~=1
    J_Musicdpd=zeros(length(yy),length(xx));
    for kk1=1:length(xx)
        for kk2=1:length(yy)
            J_Musicdpd(kk2,kk1)=-Music_dpd(Us,N_paths,B,P_bs,[xx(kk1),yy(kk2)],D1,fn(1:LenFreWindow).',lamda);
        end
    end
    J_Musicdpd=J_Musicdpd/max(max(J_Musicdpd));
    
    subplot(1,3,3);
    contour(xx,yy,J_Musicdpd);
    hold on;
    scatter3(XmirrorUE(1,1),YmirrorUE(1,1),1,100,'r','p','filled');
    scatter3(MusicP_ue(1),MusicP_ue(2),1,100,'b','+');
    
    scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
    scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
    scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
    scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');
    
    grid on;
    title({'dpd-music position spectrum';...
        ['error=' num2str(error3) 'm']});
    xlabel('x/m');
    ylabel('y/m');
    legend('Power','UE','estimated position','BS', '1#mirrorUE','2#mirrorUE','3#mirrorUE', 'Location','best');
    axis([0,30,0,30]);
end
toc