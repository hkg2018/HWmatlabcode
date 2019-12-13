function J=spotfi_FT0(Y_j,B,tau_theta,N_v,N_h,N_fft,delta_f)
a_thetax=exp(-1j*pi*(0:N_v-1)*cos(tau_theta(2))).';
a_thetay=exp(-1j*pi*(0:N_h-1)*sin(tau_theta(2))).';
a_theta=kron(a_thetay,a_thetax);
a_tau=exp(-1j*2*pi*delta_f*(0:N_fft-1)*tau_theta(1)).';
a=kron(a_tau,B*a_theta);
% J=-1/(norm(a)^2-norm(a'*Y_t)^2);
J=norm(a'*Y_j)^2;
end