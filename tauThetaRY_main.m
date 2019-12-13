% clear all;
Bw=200e6;%信号带宽
M=64;
N_fft=2048/8;
Res_aoa=1/180*pi;%多径角度间隔
Res_tau=0.8/Bw;%时延间隔
load(['..\mat\' 'ReceivedOFDMsignal_' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat']);
tic
Tscale=1e9;
Ascale=1e2;
tau_theta_e=zeros(2,N_mp,N_bs);
error=zeros(2,N_mp,N_bs);
 Us1=zeros(N_fft,N_mp,N_bs);
 Us2=zeros(M,N_mp,N_bs);
 for nbs=1:N_bs 
    RY1=zeros(N_fft);
    RY2=zeros(M);
    for q=1:Q      
        RY2=RY2+1/(Q*N_fft)*Y_t(:,:,q,nbs)*Y_t(:,:,q,nbs).';
        RY1=RY1+1/(Q*M)*Y_t(:,:,q,nbs).'*conj(Y_t(:,:,q,nbs));
    end


     [U1,tezhengzhi1,~]=svd(RY1);
     [U2,tezhengzhi2,~]=svd(RY2);
    Us1(:,:,nbs)=U1(:,1:N_mp);
    Us2(:,:,nbs)=U2(:,1:N_mp);

%     for nmp=1:N_mp
%         tau_theta0=[tau(nmp,nbs)*Tscale;theta(nmp,nbs)*Ascale];%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
%         [tau_theta_e(1,nmp,nbs),fval1,~]=fminsearch(@(tau)tau_music0(Us1(:,:,nbs),tau,N_fft,delta_f,Tscale),tau_theta0(1)) ;
%         [tau_theta_e(2,nmp,nbs),fval2,~]=fminsearch(@(theta)theta_music0(Us2(:,:,nbs),B,theta,N_v,N_h,Ascale),tau_theta0(2)) ;
%         tau_theta_e(:,nmp,nbs)=diag([1e9/Tscale 180/pi/Ascale])*tau_theta_e(:,nmp,nbs);
%         error(:,nmp,nbs)=abs(tau_theta_e(:,nmp,nbs)-[tau(nmp,nbs)*1e9;theta(nmp,nbs)*180/pi]);
%     end
 end
 toc

%% 绘制参数谱
for nbs=1:N_bs
    Rtau_e=tau(1,nbs)-Res_tau/5:Res_tau/5:tau(end,nbs)+Res_tau/5;
    Rtheta_e=(theta(2,nbs)-Res_aoa/5:Res_aoa/5:theta(end,nbs)+Res_aoa/5)/pi*180;
    J_tau=zeros(length(Rtau_e),1);
    J_theta=zeros(length(Rtheta_e),1);
    
    for kk1=1:length(Rtau_e)
        J_tau(kk1)=tau_music(Us1(:,:,nbs),Rtau_e(kk1),N_fft,delta_f);
    end
    
    for kk2=1:length(Rtheta_e)
          J_theta(kk2)=theta_music(Us2(:,:,nbs),B,Rtheta_e(kk2)/180*pi,N_v,N_h);
    end
    
    
    fig((nbs-1)*2+1)=figure;
    plot(Rtau_e,10*log10(J_tau/max(J_tau)));
    hold on;
    scatter(tau(:,nbs),ones(N_mp,1)*max(10*log10(J_tau/max(J_tau))),100,'r','p','filled');
%     scatter(tau_theta_e(1,:,nbs).'/Tscale,ones(N_mp,1)*max(10*log10(J_tau/max(J_tau))),100,'b','x');
    grid on;
    
    fig((nbs-1)*2+2)=figure;
    plot(Rtheta_e,10*log10(J_theta/max(J_theta)));
    hold on;
    scatter(theta(:,nbs)/pi*180,ones(N_mp,1)*max(10*log10(J_theta/max(J_theta))),100,'r','p','filled');
%     scatter(tau_theta_e(2,:,nbs).'/Ascale/pi*180,ones(N_mp,1)*max(10*log10(J_theta/max(J_theta))),100,'b','x');
    grid on;    
end
% savefig(fig,['..\mat\' 'tauTheta_fig' '_Resaoa' num2str(Res_aoa/pi*180) 'deg'  '_Restau' num2str(Res_tau*1e9) 'ns.fig']);