function J=tau_music0(Us,tau,N_fft,delta_f,Tscale)
a_tau=exp(-1j*2*pi*delta_f*(0:N_fft-1)*tau/Tscale).';
J=-1/(norm(a_tau)^2-norm(a_tau'*Us)^2);
end