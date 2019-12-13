function err_TOA=spotfi_Location_TOA(mea_TOA,s,P)%meaTOA是3个基站的TOA，s是基站坐标，P是目标坐标，都是按列排列
c=3e8;
delta_r1=c*(mea_TOA(2)-mea_TOA(1));
delta_r2=c*(mea_TOA(3)-mea_TOA(1));
A0 = [s(1,1)-s(1,2),s(2,1)-s(2,2);s(1,1)-s(1,3),s(2,1)-s(2,3)];
tp1=inv(A0);
k1=0.5*(delta_r1^2+(s(1,1)^2+s(2,1)^2)-(s(1,2)^2+s(2,2)^2));
k2=0.5*(delta_r2^2+(s(1,1)^2+s(2,1)^2)-(s(1,3)^2+s(2,3)^2));
m1=tp1(1,1)*k1+tp1(1,2)*k2;
n1=tp1(1,1)*delta_r1+tp1(1,2)*delta_r2;
m2=tp1(2,1)*k1+tp1(2,2)*k2;
n2=tp1(2,1)*delta_r1+tp1(2,2)*delta_r2;
A=n1^2+n2^2-1;
B=2*((m1-s(1,1))*n1+(m2-s(2,1))*n2);
C=(m1-s(1,1))^2+(m2-s(2,1))^2;
syms x
f=A*(x^2)+B*x+C;
r1=solve(f==0,x);
r0=double(r1);
u_e_TOA_1=[m1+n1*r0(1);m2+n2*r0(1)];
u_e_TOA_2=[m1+n1*r0(2);m2+n2*r0(2)];
u_e_TOA=u_e_TOA_1;
if(u_e_TOA_1(1)>35||u_e_TOA_1(1)<5||u_e_TOA_1(2)>15||u_e_TOA_1(2)<5)
u_e_TOA=u_e_TOA_2;
end
err_TOA = norm(u_e_TOA - P);
end