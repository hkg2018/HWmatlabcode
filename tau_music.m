function J=tau_music(Us,tau,N_fft,delta_f)
a_tau=exp(-1j*2*pi*delta_f*(0:N_fft-1)*tau).';
J=1/(norm(a_tau)^2-norm(a_tau'*Us)^2);
end