% clear all;
tic
%% �������
c=3e8;%����
fc=30e9;
lamda=c/fc;%�ز�����

P_ue=[15,15,4];%UEλ��
P_bs=[
    0,0,4;
    30,0,4;
    15,30,4
    ];%BSλ��
N_bs=size(P_bs,1);%BS����

N_v=8;%��������
N_h=8;%��������

d=lamda/2;%��Ԫ���
M=N_h*N_v;%������

D0=zeros(3,M);%����Ԫ�����λ������(BS������״)
for nh=1:N_h
    for nv=1:N_v
        D0(:,(nh-1)*N_v+nv)=[(nv-1)*d;(nh-1)*d;0];
    end
end

delta_f=120e3*8;%OFDMƵ�ʼ��
N_fft=2048/8;%FFT����
Bw=200e6;%�źŴ���

% load(['..\mat\' 'B0']);
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);��ͨ�����������뱾
gama0=[(pi/4) (3*pi/4) (-pi/2)];%UE����BS�ľ�������Ƕ�
gama_H0=gama0-[(0),(pi/4),(0)];%�����������е���ת�Ƕ�=��������Ƕ�-�������Ƕ�
Res_aoa=0.5/180*pi;%�ྶ�Ƕȼ��
% range_aoaH=[-fix((N_mp-1)/2):-1 1:ceil((N_mp-1)/2)]*Res_aoa;%����վ�Ķྶ����Ƕ�
% range_aoaV=[-fix((N_mp-1)/2):-1 1:ceil((N_mp-1)/2)]*Res_aoa;%�߶ȽǷ�Χ
range_aoaH0=-pi/3:Res_aoa:pi/3;%����վ�ķ�λ�Ƿ�Χ
range_aoaV0=-pi/6:Res_aoa:pi/6;%�߶ȽǷ�Χ

SNR=20;%�����dB��%ĿǰSNR������symbol�����������������֮��
sigma2=10^(-174/10-3)*Bw;%�������ʣ�W
Q=20;%������

%����������
CFO=0;%�ز�Ƶ��ƫ��1.5e3
STO=0;%������ʱƫ��
SFO=0;%����Ƶ��Ƶ��µĲ�������ƫ��0.6e-9/10e-3

T=1/delta_f;%�۲�ʱ��s,����Ƶ�ʷֱ���
fs=N_fft*delta_f;%������������ʱ��ֱ���
Ts=1/fs;%�źŲ�������
N_ofdm=fix(Bw/delta_f);%OFDM���ز���
fn0=(1:N_ofdm)*delta_f;%�����ź�Ƶ��
tt=(0:N_fft-1)*Ts;%BSʱ���ᣬ��λ��


%% ���ɽ����ź�

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

Y_t=zeros(M,N_fft,Q,N_bs);%ͨ�����Ƶ������
parfor q=1:Q
    
    alpha=alpha0;%ֱ�ﾶ�ŵ������һ��
    tau=tau0/1e-9+sqrt(diag((P_ue-P_bs)*(P_ue-P_bs).')).'/c;%�ŵ��弤��Ӧ������ֱ�ﲨʱ��
    N_mp=N_mp0;
    range_aoaH=range_aoaH0;
    gama_H=gama_H0;
    gama=gama0;
    P_bs0=P_bs;
    fn=fn0;
    Xk=sqrt(10^(SNR/10)*sigma2)*ones(1,N_ofdm);%�ŵ���������;   exp(1j*2*pi*rand(1,N_ofdm));%�����Ƶ����
    
    for nbs=1:N_bs
        
        phi=zeros(N_mp(nbs)-1,1);%�ྶ�߶Ƚ�
        F_alpha=ones(N_mp(nbs),1);%�ྶϵ������˥������,��ʼ��Ϊ��˥��
        
        theta=[0;range_aoaH(randi(length(range_aoaH),1,N_mp(nbs)-1)).']+gama(nbs)-gama_H(nbs);
        %    pp_aoaV=randperm(fix(pi/3/Res_aoa));
        %    phi=-pi/6+pp_aoaV(1:N_mp-1)*Res_aoa;
        Kvec=-2*pi/lamda*[(P_ue-P_bs0(nbs,:))/norm(P_ue-P_bs0(nbs,:));(cos(phi).*cos(gama_H(nbs)+theta(2:end,1))) (cos(phi).*sin(gama_H(nbs)+theta(2:end,1))) sin(phi)];
        
        D=[
            cos(gama_H(nbs)),sin(-gama_H(nbs)),0;
            -sin(-gama_H(nbs)),cos(gama_H(nbs)),0;
            0,0,1
            ]*D0;%��ת��������
        
        A=exp(1j*(Kvec*D).');
        
        %F_alpha(:,nbs)=exp((tau(1,nbs)-tau(:,nbs))/tau_max);
        
        
        %�����ض������½����ź��е�ȷ������

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