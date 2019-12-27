% clear all;
tic
%% �������
c=3e8;%����
lamda=c/30e9;%�ز�����

P_bs=[
    0,0,4;
    30,0,4;
    15,30,4
    ];%BSλ��
N_bs=size(P_bs,1);%BS����

N_v=8;%��������
N_h=8;%��������
d=lamda/2;%��Ԫ���
M=N_h*N_v;%������
gama_H0=[pi/4,3*pi/4,-pi/2];%�����������е���ת�Ƕ�=��������Ƕ�-�������Ƕ�

D0=zeros(3,M);%����Ԫ�����λ������(BS������״)
D1=zeros(3,M,N_bs);
for nh=1:N_h
    for nv=1:N_v
        D0(:,(nh-1)*N_v+nv)=[(nv-1)*d;(nh-1)*d;0];
    end
end
for nbs=1:N_bs
    D1(:,:,nbs)=[
        cos(gama_H0(nbs)),sin(-gama_H0(nbs)),0;
        -sin(-gama_H0(nbs)),cos(gama_H0(nbs)),0;
        0,0,1
        ]*D0;%��ת��������
end




delta_f=120e3*8;%OFDMƵ�ʼ��
N=207;%FFT����

M_beam=size(B,1);
fn=(1:N)*delta_f;%�����ź�Ƶ��
B=DFTcodebook_spotfimuci_1204(N_v,N_h);%eye(M);��ͨ�����������뱾
LenFreWindow=fix(N*0.25);
Qf=N-LenFreWindow+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ydpd_t=zeros(N*M_beam,N_bs);%ͨ�����Ƶ������
P_ue=[15,15,4];%UEλ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Us=zeros(LenFreWindow*M_beam,20,N_bs);%ͨ�����Ƶ������
tez=zeros(:,N_bs);
N_paths=zeros(1,N_bs);
for nbs=1:N_bs
Y_ph=zeros(LenFreWindow*M_beam,Qf);
Y0=Ydpd_t(:,:,nbs);

for qf=1:Qf
    Y_ph(:,qf)=Y0((qf-1)*M_beam+(1:LenFreWindow*M_beam),:);
end
[U,tez0,~]=svd(Y_ph);
tez(:,nbs)=diag(tez0);
[IDX,~]=kmeans(tez(:,nbs),2);
N_paths(nbs)=length(find(IDX==IDX(1)));
Us(:,:,nbs)=U(:,1:N_paths(nbs));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%λ�ù���%%%%%%%%%%%%%%%%%%%%%%%%%%
%             P_ue20=P_ue(1:2);%+[Res_tau/4*Tscale;Res_aoa/4*Ascale]*rand(1,1);

[MusicP_ue,~,~]=fminsearch(@(P_ue_xy)Music_dpd(Us,N_paths,B,P_bs,P_ue_xy,P_ue_z,D1,fn(1:LenFreWindow).',lamda,h_re,tau_re),P_ue(1:2));
% [MF1P_ue,~,~]=fminsearch(@(P_ue2)MF1dpd(Ydpd_t,B,P_bs,P_ue2,D1,fn.',lamda),P_ue(1:2));
[MF0P_ue,~,~]=fminsearch(@(P_ue_xy)MF0dpd(Ydpd_t,B,P_bs,P_ue_xy,P_ue_z,D1,fn,lamda,h_re,tau_re),MF1P_ue);

error1=norm([MF0P_ue,4]-P_ue);
% error2=norm([MF1P_ue,4]-P_ue);
error3=norm([MusicP_ue,4]-P_ue);
toc

%% ����λ����
tic
xx=1:1:29;
yy=1:1:29;
% J_dpdMF1=zeros(length(yy),length(xx));
J_dpdMF0=zeros(length(yy),length(xx));

for kk1=1:length(xx)
    for kk2=1:length(yy)
        J_dpdMF0(kk2,kk1)=-MFdpd(Ydpd_t,B,P_bs,[xx(kk1),yy(kk2)],D1,fn.',lamda);
%         J_dpdMF1(kk2,kk1)=-FastHR_dpd(Ydpd_t,B,P_bs,[xx(kk1),yy(kk2)],D1,fn.',lamda);
    end
end
J_dpdMF0=J_dpdMF0/max(max(J_dpdMF0));
% J_dpdMF1=J_dpdMF1/max(max(J_dpdMF1));


% [XmirrorUE,YmirrorUE]=pol2cart(theta+gama_H0,tau0*c);
% XmirrorUE=XmirrorUE+P_bs(:,1).';
% YmirrorUE=YmirrorUE+P_bs(:,2).';

figure('NumberTitle', 'off', 'Name', ...
    [ 'N_{paths}=' num2str(N_paths) ',Tau-Theta_{res}[' num2str(min_tau) 'ns,' num2str(min_aoa) 'deg]' 'CFO' num2str(CFO0) 'Hz']);
subplot(1,3,1);
contour(xx,yy,J_dpdMF0);
hold on;
scatter3(P_ue(1),P_ue(2),1,100,'r','p','filled');
scatter3(MF0P_ue(1),MF0P_ue(2),1,100,'b','+');

scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
% scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
% scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
% scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');

grid on;
title({'dpd-MF position spectrum';...
    ['error=' num2str(error1) 'm']});
xlabel('x/m');
ylabel('y/m');

axis([0,30,0,30]);

subplot(1,3,2);
contour(xx,yy,J_dpdMF1);
hold on;
scatter3(P_ue(1),P_ue(2),1,100,'r','p','filled');
scatter3(MF1P_ue(1),MF1P_ue(2),1,100,'b','+');

scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
% scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
% scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
% scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');
grid on;
title({'improved dpd-MF position spectrum';...
    ['error=' num2str(error2) 'm']});
xlabel('x/m');
ylabel('y/m');
axis([0,30,0,30]);




    J_Musicdpd=zeros(length(yy),length(xx));
    for kk1=1:length(xx)
        for kk2=1:length(yy)
            J_Musicdpd(kk2,kk1)=-Music_dpd(Us,N_paths,B,P_bs,[xx(kk1),yy(kk2)],D1,fn(1:LenFreWindow).',lamda);
        end
    end
    J_Musicdpd=J_Musicdpd/max(max(J_Musicdpd));
    
    subplot(1,3,3);
    contour(xx,yy,J_Musicdpd);
    hold on;
    scatter3(P_ue(1),P_ue(2),1,100,'r','p','filled');
    scatter3(MusicP_ue(1),MusicP_ue(2),1,100,'b','+');
    
    scatter3(P_bs(:,1),P_bs(:,2),ones(N_bs,1),100,'b','s','filled');
%     scatter3(XmirrorUE(2:N_paths(1),1),YmirrorUE(2:N_paths(1),1),ones(N_paths(1)-1,1),50,'r','p');
%     scatter3(XmirrorUE(2:N_paths(2),2),YmirrorUE(2:N_paths(2),2),ones(N_paths(2)-1,1),50,'b','p');
%     scatter3(XmirrorUE(2:N_paths(3),3),YmirrorUE(2:N_paths(3),3),ones(N_paths(3)-1,1),50,'g','p');
    
    grid on;
    title({'dpd-music position spectrum';...
        ['error=' num2str(error3) 'm']});
    xlabel('x/m');
    ylabel('y/m');
    legend('Power','UE','estimated position','BS', 'Location','best');
    axis([0,30,0,30]);

toc