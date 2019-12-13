tic
Tscale=1e9;
Ascale=1e2;
tau_theta_e=zeros(2,N_paths,N_bs);
error=zeros(2,N_paths,N_bs);
Music_flag=0;
for nbs=1:N_bs

    for nmp=1:N_paths
        tau_theta0=[tau0(nmp,nbs)*Tscale;theta(nmp,nbs)*Ascale];%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
        
        if Music_flag
            [tau_theta_e(:,nmp,nbs),~,~]=fminsearch(@(tau_theta)spotfi_bs(Us(:,:,nbs),B,tau_theta,N_v,N_h,LenFreWindow,delta_f,Tscale,Ascale),tau_theta0);
        else
            [tau_theta_e(:,nmp,nbs),~,~]=fminsearch(@(tau_theta)spotfi_FT(Ydpd_t(:,:,nbs),B,tau_theta,N_v,N_h,N_fft,delta_f,Tscale,Ascale),tau_theta0);
        end
        
        tau_theta_e(:,nmp,nbs)=diag([1e9/Tscale 180/pi/Ascale])*tau_theta_e(:,nmp,nbs);
    end
    [~,I]=sort(tau_theta_e(1,:,nbs));
    tau_theta_e(:,:,nbs)=tau_theta_e(:,I,nbs);
    error(:,:,nbs)=abs(tau_theta_e(:,:,nbs)-[tau0(:,nbs).'*1e9;theta(:,nbs).'*180/pi]);    
end
toc
% save(['..\mat\' 'spotfi_Us' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat'],'Us');
% save(['..\mat\' 'spotfi_Us' num2str(Res_aoa/pi*180) 'deg_' 'IEEE 802.15.4aChannel model'  'M' num2str(M*N_fft) 'Nfft.mat'],'Us');

%% »æÖÆ²ÎÊýÆ×
tic
figure;
if ~Music_flag
    for nbs=1:N_bs
        Rtau_e=(min(tau0(:,nbs))-2*Res_tau:Res_tau/10:max(tau0(:,nbs))+Res_tau)*1e9;
        Rtheta_e=(min(theta(:,nbs))-Res_aoa:Res_aoa/10:max(theta(:,nbs))+Res_aoa)/pi*180;
        J_spotfiFT=zeros(length(Rtheta_e),length(Rtau_e));
        for kk1=1:length(Rtau_e)
            for kk2=1:length(Rtheta_e)
                J_spotfiFT(kk2,kk1)=spotfi_FT0(Ydpd_t(:,:,nbs),B,[Rtau_e(kk1)/1e9;Rtheta_e(kk2)/180*pi],N_v,N_h,N_fft,delta_f);
            end
        end
        J_spotfiFT=10*log10(J_spotfiFT/max(max(J_spotfiFT)));
        subplot(1,3,nbs);
        contour(Rtau_e,Rtheta_e,J_spotfiFT);
        hold on;
        scatter3(tau0(1,nbs)*1e9,theta(1,nbs)/pi*180,0,100,'r','p','filled');
        scatter3(tau_theta_e(1,:,nbs).',tau_theta_e(2,:,nbs).',zeros(N_paths,1),100,'r','x');
        scatter3(tau0(2:end,nbs)*1e9,theta(2:end,nbs)/pi*180,zeros(N_paths-1,1),100,'b','p','filled');
        %     scatter3(tau_theta_e(1,2:end,nbs).'/1e9,tau_theta_e(2,2:end,nbs).',zeros(N_paths-1,1),100,'b','x');
        grid on;
        title({ [num2str(nbs) '#spotfiSpectrum'];[ 'N_{paths}=' num2str(N_paths) ',Tau-Theta_{res}[' num2str(Res_tau*1e9) 'ns,' num2str(Res_aoa/pi*180) 'deg],'];...
            ['tauError(ns)=' num2str(error(1,1,nbs))];['thetaError(deg)=' num2str(error(2,1,nbs))]});
        xlabel('tau/ns');
        ylabel('theta/deg');
        if nbs==1
            legend('Power','true-tauTheta','estimated-tauTheta','true-tauThetaMP', 'Location','best');
        end
    end
else
    for nbs=1:N_bs
        Rtau_e=(min(tau0(:,nbs))-2*Res_tau:Res_tau/5:max(tau0(:,nbs))+Res_tau)*1e9;
        Rtheta_e=(min(theta(:,nbs))-Res_aoa:Res_aoa/5:max(theta(:,nbs))+Res_aoa)/pi*180;
        J_spotfiMusic=zeros(length(Rtheta_e),length(Rtau_e));
        for kk1=1:length(Rtau_e)
            for kk2=1:length(Rtheta_e)
                J_spotfiMusic(kk2,kk1)=spotfi_bs0(Us(:,:,nbs),B,[Rtau_e(kk1)/1e9;Rtheta_e(kk2)/180*pi],N_v,N_h,LenFreWindow,delta_f);
            end
        end
        J_spotfiMusic=10*log10(J_spotfiMusic/max(max(J_spotfiMusic)));
        subplot(1,3,nbs);
        contour(Rtau_e,Rtheta_e,J_spotfiMusic);
        hold on;
        scatter3(tau0(1,nbs)*1e9,theta(1,nbs)/pi*180,0,100,'r','p','filled');
        scatter3(tau_theta_e(1,1,nbs).',tau_theta_e(2,1,nbs).',0,100,'r','x');
        scatter3(tau0(2:end,nbs)*1e9,theta(2:end,nbs)/pi*180,zeros(N_paths-1,1),100,'b','p','filled');
        %     scatter3(tau_theta_e(1,2:end,nbs).'/1e9,tau_theta_e(2,2:end,nbs).',zeros(N_paths-1,1),100,'b','x');
        grid on;
        title({ [num2str(nbs) '#spotfiSpectrum'];[ 'N_{paths}=' num2str(N_paths) ',Tau-Theta_{res}[' num2str(Res_tau*1e9) 'ns,' num2str(Res_aoa/pi*180) 'deg],'];...
            ['tauError(ns)=' num2str(error(1,1,nbs))];['thetaError(deg)=' num2str(error(2,1,nbs))]});
        xlabel('tau/ns');
        ylabel('theta/deg');
        if nbs==1
            legend('Power','true-tauTheta','estimated-tauTheta','true-tauThetaMP', 'Location','best');
        end
    end
end
toc
%  savefig(fig,['..\mat\' 'spotfi_fig' '_Resaoa' num2str(Res_aoa/pi*180) 'deg'  '_Restau' num2str(Res_tau*1e9) 'ns.fig']);