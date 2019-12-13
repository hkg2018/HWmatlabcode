% clear all;
tic
% Y_t0=permute(Y_t,[1,3,2,4]);
Tscale=1;
Ascale=1;
tau_theta_e=zeros(2,N_paths,N_bs);
error=zeros(2,N_paths,N_bs);

 Us1=zeros(N_fft,N_paths,N_bs);
 Us2=zeros(M,N_paths,N_bs);
 for nbs=1:N_bs 
    RY1=zeros(N_fft);
    RY2=zeros(M);
    for q=1:Q      
        RY2=RY2+1/(Q*N_fft)*(Y_t0(:,:,q,nbs)*Y_t0(:,:,q,nbs)');
        Y0=Y_t0(:,:,q,nbs).';
        RY1=RY1+1/(Q*M)*(Y0*Y0');
    end
%     U2=zeros(M);
%     U1=zeros(N_fft);

% %     [U2,qiyizhi,~]=svd(Y1(:,:,nbs));
% %     [U1,qiyizhi1,~]=svd(Y1(:,:,nbs).');
    [U1,tezhengzhi1,~]=svd(RY1);
    [U2,tezhengzhi2,~]=svd(RY2);
    Us1(:,:,nbs)=U1(:,1:N_paths);
    Us2(:,:,nbs)=U2(:,1:N_paths);

    for nmp=1:N_paths
        tau_theta0=[tau0(nmp,nbs)*Tscale;theta(nmp,nbs)*Ascale];%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
        [tau_theta_e(1,nmp,nbs),fval1,~]=fminsearch(@(tau)tau_music0(Us1(:,:,nbs),tau,N_fft,delta_f,Tscale),tau_theta0(1)) ;
        [tau_theta_e(2,nmp,nbs),fval2,~]=fminsearch(@(theta)theta_music0(Us2(:,:,nbs),B,theta,N_v,N_h,Ascale),tau_theta0(2)) ;
        tau_theta_e(:,nmp,nbs)=diag([1e9/Tscale 180/pi/Ascale])*tau_theta_e(:,nmp,nbs);             
    end
    tau_theta_e(1,:,nbs)=sort(tau_theta_e(1,:,nbs));    
    error(:,:,nbs)=abs(tau_theta_e(:,:,nbs)-[tau0(:,nbs).'*1e9;theta(:,nbs).'*180/pi]);
 end
 toc

%% »æÖÆ²ÎÊýÆ×
figure;
for nbs=1:N_bs
    Rtau_e=min(tau0(:,nbs))-2*Res_tau:0.1*1e-9:max(tau0(:,nbs))+Res_tau;
    Rtheta_e=(min(theta(:,nbs))-Res_aoa:0.01/180*pi:max(theta(:,nbs))+Res_aoa)/pi*180;
    J_tau=zeros(length(Rtau_e),1);
    J_theta=zeros(length(Rtheta_e),1);
    
    for kk1=1:length(Rtau_e)
        J_tau(kk1)=tau_music(Us1(:,:,nbs),Rtau_e(kk1),N_fft,delta_f);
    end
    
    for kk2=1:length(Rtheta_e)
          J_theta(kk2)=theta_music(Us2(:,:,nbs),B,Rtheta_e(kk2)/180*pi,N_v,N_h);
    end
    
    J_tau=10*log10(J_tau/max(J_tau));
    subplot(2,3,nbs);
    plot(Rtau_e*1e9,J_tau);
    hold on;
    scatter(tau0(1,nbs)*1e9,0,100,'r','p','filled');
    scatter(tau_theta_e(1,1,nbs).',0,100,'r','x');
    scatter(tau0(2:end,nbs)*1e9,zeros(N_paths-1,1),100,'b','p','filled');
%     scatter(tau_theta_e(1,2:end,nbs).'/1e9,zeros(N_paths-1,1),100,'b','x');
    grid on;
    title({[num2str(nbs) '#tauMusicSpectrum'];[ 'N_{paths}=' num2str(N_paths) ',Tau_{res}=' num2str(Res_tau*1e9) 'ns,'];...
       [ 'Error(ns)=' num2str(error(1,1,nbs))]});
    xlabel('tau/ns');
    ylabel('power/dB');
    
    if nbs==1
    legend('spectrum','true-tau','estimated-tau','true-tauMP', 'Location','best');
    end
    
    J_theta=10*log10(J_theta/max(J_theta));
    subplot(2,3,nbs+3);
    plot(Rtheta_e,J_theta);
    hold on;
    scatter(theta(1,nbs)/pi*180,0,100,'r','p','filled');
    scatter(tau_theta_e(2,1,nbs).',0,100,'r','x');
    scatter(theta(2:end,nbs)/pi*180,zeros(N_paths-1,1),100,'b','p','filled');
%     scatter(tau_theta_e(2,2:end,nbs).',zeros(N_paths-1,1),100,'b','x');
    grid on;    
        title({[num2str(nbs) '#thetaMusicSpectrum'];[ 'N_{paths}=' num2str(N_paths) ',Theta_{res}=' num2str(Res_aoa/pi*180) 'deg,'];...
        ['Error(deg)=' num2str(error(2,1,nbs))]});
    xlabel('theta/deg');
    ylabel('power/dB'); 
    if nbs==1
    legend('spectrum','true-theta','estimated-theta','true-thetaMP', 'Location','best');   
    end
end

% savefig(fig,['..\mat\' 'tauTheta_fig' '_Resaoa' num2str(Res_aoa/pi*180) 'deg'  '_Restau' num2str(Res_tau*1e9) 'ns.fig']);