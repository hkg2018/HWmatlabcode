function [tau_e,theta_e]=tauTheta_MUSICestimator(Y_t,B,N_v,N_h,N_fft,delta_f,tau0,theta,N_paths,gama_H)
Y_t0=permute(Y_t,[1,3,2,4]);
M=size(Y_t0,1);
Q=size(Y_t0,3);
Tscale=1e9;
Ascale=1e2;
N_bs=size(Y_t0,4);
tau_theta_e=zeros(2,N_bs);
tau_theta_e0=zeros(2,N_paths);

 for nbs=1:N_bs 
     RY1=zeros(N_fft);
     RY2=zeros(M);
    for q=1:Q      
        RY2=RY2+1/(Q*N_fft)*(Y_t0(:,:,q,nbs)*Y_t0(:,:,q,nbs)');
        Y0=Y_t0(:,:,q,nbs).';
        RY1=RY1+1/(Q*M)*(Y0*Y0');
    end

    [U1,~,~]=svd(RY1);
    [U2,~,~]=svd(RY2);
    Us1=U1(:,1:N_paths);
    Us2=U2(:,1:N_paths);
    fval=zeros(N_paths,1);
    for nmp=1:N_paths
        tau_theta0=[tau0(nmp,nbs)*Tscale;theta(nmp,nbs)*Ascale];%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);
        [tau_theta_e0(1,nmp),~,~]=fminsearch(@(tau)tau_music0(Us1,tau,N_fft,delta_f,Tscale),tau_theta0(1)) ;
        [tau_theta_e0(2,nmp),fval(nmp),~]=fminsearch(@(theta)theta_music0(Us2,B,theta,N_v,N_h,Ascale),tau_theta0(2)) ;
        tau_theta_e0(:,nmp)=diag([1/Tscale 180/pi/Ascale])*tau_theta_e0(:,nmp);             
    end
    tau_theta_e(1,nbs)=min(tau_theta_e0(1,:));
     [~,I]=min(fval);
    tau_theta_e(2,nbs)=tau_theta_e0(2,I)+gama_H(nbs)/pi*180;   
 end
 tau_e=tau_theta_e(1,:);
theta_e=tau_theta_e(2,:);
end







