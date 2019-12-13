function J=theta_music(Us,B,theta,N_v,N_h)
a_thetax=exp(-1j*pi*(0:N_v-1)*cos(theta)).';
a_thetay=exp(-1j*pi*(0:N_h-1)*sin(theta)).';
a_theta=B*kron(a_thetay,a_thetax);
% J=1/(norm(a_theta)^2-norm(a_theta'*Us)^2);
J=1/(norm(a_theta)^2-norm(a_theta'*Us)^2);
end