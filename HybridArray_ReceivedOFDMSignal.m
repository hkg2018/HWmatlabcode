% clear all;
tic
%% 仿真参数
c=3e8;%光速
fc=30e9;
lamda=c/fc;%载波波长

P_ue=[15,15,4];%UE位置
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

D0=zeros(3,M);%各阵元的相对位置坐标(BS阵列形状)
for nh=1:N_h
    for nv=1:N_v
        D0(:,(nh-1)*N_v+nv)=[(nv-1)*d;(nh-1)*d;0];
    end
end

delta_f=120e3*8;%OFDM频率间隔
N_fft=2048/8;%FFT点数
Bw=200e6;%信号带宽

% load(['..\mat\' 'B0']);
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);单通道正交波束码本
gama0=[(pi/4) (3*pi/4) (-pi/2)];%UE到各BS的绝对入射角度
gama_H0=gama0-[(0),(pi/4),(0)];%各个天线阵列的旋转角度=绝对入射角度-相对入射角度
Res_aoa=0.5/180*pi;%多径角度间隔
% range_aoaH=[-fix((N_mp-1)/2):-1 1:ceil((N_mp-1)/2)]*Res_aoa;%各基站的多径入射角度
% range_aoaV=[-fix((N_mp-1)/2):-1 1:ceil((N_mp-1)/2)]*Res_aoa;%高度角范围
range_aoaH0=-pi/3:Res_aoa:pi/3;%各基站的方位角范围
range_aoaV0=-pi/6:Res_aoa:pi/6;%高度角范围

SNR=20;%信噪比dB，%目前SNR定义是symbol功率与测量噪声功率之比
sigma2=10^(-174/10-3)*Bw;%噪声功率，W
Q=20;%样本数

%非理想因素
CFO=0;%载波频率偏差1.5e3
STO=0;%采样定时偏差
SFO=0;%采样频率频差导致的采样周期偏差0.6e-9/10e-3

T=1/delta_f;%观测时长s,决定频率分辨率
fs=N_fft*delta_f;%基带带宽，决定时间分辨率
Ts=1/fs;%信号采样周期
N_ofdm=fix(Bw/delta_f);%OFDM子载波数
fn0=(1:N_ofdm)*delta_f;%基带信号频点
tt=(0:N_fft-1)*Ts;%BS时间轴，单位秒


%% 生成接收信号

%%%%%%%%%%%%%%%%%%%%%%%%%%IEEE 802.15.4aChannel model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_channels = N_bs; % number of channel impulse responses to generate
rng(12);% initialize state of function for repeatability

cm_num = 3; % channel model number from 1 to 8
% get channel model params based on this channel model number
[Lam,lambda,Lmean,lambda_mode,lambda_1,lambda_2,beta,Gam,gamma_0,Kgamma, ...
    sigma_cluster,nlos,gamma_rise,gamma_1,chi,m0,Km,sigma_m0,sigma_Km, ...
    sfading_mode,m0_sp,std_shdw,kappa,~,~] = uwb_sv_params_15_4a( cm_num );

% get a bunch of realizations (impulse responses)
[alpha0,tau0,t0,N_mp0] = uwb_sv_model_ct_15_4a(Lam,lambda,Lmean,lambda_mode,lambda_1, ...
    lambda_2,beta,Gam,gamma_0,Kgamma,sigma_cluster,nlos,gamma_rise,gamma_1, ...
    chi,m0,Km,sigma_m0,sigma_Km,sfading_mode,m0_sp,std_shdw,num_channels,Ts*1e9);

alpha0=alpha0./alpha0(1,:).*exp(2j*pi*rand(size(alpha0)));
%%%%%%%%%%%%%%%%%%%%%%%%%%IEEE 802.15.4aChannel model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_t=zeros(M,N_fft,Q,N_bs);%通道输出频域数据
parfor q=1:Q
    
    alpha=alpha0;%直达径信道增益归一化
    tau=tau0/1e-9+sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).')).'/c;%信道冲激响应中引入直达波时延
    N_mp=N_mp0;
    range_aoaH=range_aoaH0;
    gama_H=gama_H0;
    gama=gama0;
    P_bs0=P_bs;
    fn=fn0;
    Xk=sqrt(10^(SNR/10)*sigma2)*ones(1,N_ofdm);%信道测试序列;   exp(1j*2*pi*rand(1,N_ofdm));%随机导频序列
    
    for nbs=1:N_bs
        
        phi=zeros(N_mp(nbs)-1,1);%多径高度角
        F_alpha=ones(N_mp(nbs),1);%多径系数方差衰减因子,初始化为不衰减
        
        theta=[0;range_aoaH(randi(length(range_aoaH),1,N_mp(nbs)-1)).']+gama(nbs)-gama_H(nbs);
        %    pp_aoaV=randperm(fix(pi/3/Res_aoa));
        %    phi=-pi/6+pp_aoaV(1:N_mp-1)*Res_aoa;
        Kvec=-2*pi/lamda*[(P_ue-P_bs0(nbs,:))/norm(P_ue-P_bs0(nbs,:));(cos(phi).*cos(gama_H(nbs)+theta(2:end,1))) (cos(phi).*sin(gama_H(nbs)+theta(2:end,1))) sin(phi)];
        
        D=[
            cos(gama_H(nbs)),sin(-gama_H(nbs)),0;
            -sin(-gama_H(nbs)),cos(gama_H(nbs)),0;
            0,0,1
            ]*D0;%旋转天线阵列
        
        A=exp(1j*(Kvec*D).');
        
        %F_alpha(:,nbs)=exp((tau(1,nbs)-tau(:,nbs))/tau_max);
        
        
        %生成特定场景下接收信号中的确定部分

        m_t0=zeros(M,N_fft);
        for nmp=1:N_mp(nbs)
            s_t=zeros(1,N_fft);
            for k=1:N_ofdm
                s_t=s_t+Xk(k)*exp(1j*2*pi*fn(k)*(tt*(1+SFO)-STO-tau(nmp,nbs)));
            end
            m_t0=m_t0+alpha(nmp,nbs)*A(:,nmp)*s_t;
        end
        m_t=B*m_t0*diag(exp(-1j*2*pi*CFO*tt));
        Y_t(:,:,q,nbs)=fft(m_t,[],2)+wgn(M,N_fft,sigma2,'linear','complex');
    end
end
toc
% save(['..\mat\' 'ReceivedOFDMsignal_' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat']);
 save(['..\mat\' 'ReceivedOFDMsignal_' num2str(Res_aoa/pi*180) 'deg_' 'IEEE 802.15.4aChannel model'  'M' num2str(M*N_fft) 'Nfft.mat']);