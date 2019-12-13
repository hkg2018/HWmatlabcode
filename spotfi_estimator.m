function [tau_e,theta_e]=spotfi_estimator(Us,B,N_v,N_h,delta_f,tau,theta,N_paths,gama_H)
Tscale=1e9;
Ascale=1e2;
M=N_v*N_h;
LenFreWindow=size(Us,1)/M;
N_bs=size(Us,3);

tau_theta_e0=zeros(2,N_paths);
tau_theta_e=zeros(2,N_bs);

for nbs=1:N_bs

    for nmp=1:N_paths
        tau_theta0=[tau(nmp,nbs)*Tscale;theta(nmp,nbs)*Ascale];%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
        [tau_theta_e0(:,nmp),~,~]=fminsearch(@(tau_theta)spotfi_bs(Us(:,:,nbs),B,tau_theta,N_v,N_h,LenFreWindow,delta_f,Tscale,Ascale),tau_theta0);
    end
    [~,I]=min(tau_theta_e0(1,:));
    tau_theta_e(:,nbs)=diag([1/Tscale 180/pi/Ascale])*tau_theta_e0(:,I)+[0;gama_H(nbs)/pi*180];
end
tau_e=tau_theta_e(1,:);
theta_e=tau_theta_e(2,:);
end


