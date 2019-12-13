% clear all;
% Bw=200e6;%信号带宽
% M=64;
% N_fft=2048/8;
% Res_aoa=1/180*pi;%多径角度间隔
% Res_tau=0.8/Bw;%时延间隔

% load(['..\mat\' 'dpdReceivedOFDMsignal_' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat']);
% load(['..\mat\' 'ReceivedOFDMsignal_' num2str(Res_aoa/pi*180) 'deg_' 'IEEE 802.15.4aChannel model'  'M' num2str(M*N_fft) 'Nfft.mat']);
tic

% load(['..\mat\' 'MVDR_dpd_invRk' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat']);
% load(['..\mat\' 'MVDR_dpd_invRk' num2str(Res_aoa/pi*180) 'deg_' 'IEEE 802.15.4aChannel model'  'M' num2str(M*N_fft) 'Nfft.mat']);

% invRk=zeros(N_bs*M,N_bs*M,N_ofdm);
% for k=1:N_ofdm
%     yk=zeros(N_bs*M,Q);
%     for q=1:Q
%         Y0=Ydpd_t(:,:,k,q);%这里是转置不是共轭转置要注意
%         yk(:,q)= Y0(:);%得到MN*MN的自相关矩阵
%     end
%     invRk(:,:,k)=pinv(1/Q*(yk*yk'));
% end

Yq=zeros(N_ofdm*N_bs*M,Q);
for k=1:N_ofdm
    yk=zeros(N_bs*M,Q);
    for q=1:Q
        Y0=Ydpd_t(:,:,k,q);%这里是转置不是共轭转置要注意
        yk(:,q)= Y0(:);%得到MN*MN的自相关矩阵
    end
    Yq((k-1)*N_bs*M+1:k*N_bs*M,:)=yk;
end


P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
[P_ue_e,fval,~]=fminsearch(@(P_ue2)MVDR_dpd(Yq,X,B,P_bs,P_ue2,D,fn,lamda),P_ue20);
error=norm([P_ue_e,4]-P_ue);
toc
% save(['..\mat\' 'MVDR_dpd_invRk' num2str(Res_aoa/pi*180) 'deg_'  num2str(Res_tau*1e9) 'ns_' 'M' num2str(M*N_fft) 'Nfft.mat'],'invRk');
% save(['..\mat\' 'MVDR_dpd_invRk' num2str(Res_aoa/pi*180) 'deg_' 'IEEE 802.15.4aChannel model'  'M' num2str(M*N_fft) 'Nfft.mat'],'invRk');

%% 绘制位置谱
tic
xx=1:1:29;
yy=1:1:29;
    J_MVDRdpd=zeros(length(yy),length(xx));
    for kk1=1:length(xx)
        for kk2=1:length(yy)
            J_MVDRdpd(kk2,kk1)=-MVDR_dpd(Yq,X,B,P_bs,[xx(kk1),yy(kk2)],D,fn,lamda);
        end
    end
    J_MVDRdpd=J_MVDRdpd/max(max(J_MVDRdpd));
    fig=figure;
    contour(xx,yy,J_MVDRdpd);
    hold on;
    scatter3(P_ue(1),P_ue(2),1,100,'r','p','filled');
    scatter3(P_ue_e(1),P_ue_e(2),1,100,'b','x');
    grid on;
    
toc
%  savefig(fig,['..\mat\' 'spotfi_fig' '_Resaoa' num2str(Res_aoa/pi*180) 'deg'  '_Restau' num2str(Res_tau*1e9) 'ns.fig']);