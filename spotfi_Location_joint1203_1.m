function Poserr=spotfi_Location_joint1203_1(canshu,P_bs,P_ue)%
%% 参数设置 
M = 3;                                    % 基站个数
C=3e8;
s=P_bs(:,1:2)';%3个基站的坐标
P=P_ue(1:2)';%目标真实坐标 
%% 联合估计联合解算+分开估计联合解算+只测向解算
      for m = 1:M
            b_joint_joint(:,m) = [cosd(canshu(1,m));sind(canshu(1,m))];
            gm_joint_joint(:,m) = [sind(canshu(1,m));-cosd(canshu(1,m))];
            b_jointDFT_joint(:,m) = [cosd(canshu(5,m));sind(canshu(5,m))];
            gm_jointDFT_joint(:,m) = [sind(canshu(5,m));-cosd(canshu(5,m))];
            b_sperate_joint(:,m) = [cosd(canshu(3,m));sind(canshu(3,m))];
            gm_sperate_joint(:,m) = [sind(canshu(3,m));-cosd(canshu(3,m))];
      end
        
        h1_joint_joint = [];h2_joint_joint = [];
        G1_joint_joint = [];G2_joint_joint = [];
        h1_jointDFT_joint = [];h2_jointDFT_joint = [];
        G1_jointDFT_joint = [];G2_jointDFT_joint = [];
        h1_sperate_joint = [];h2_sperate_joint = [];
        G1_sperate_joint = [];G2_sperate_joint = [];
        for m = 2:M
            h1_joint_joint = [h1_joint_joint;(b_joint_joint(:,m)-b_joint_joint(:,1))'*(s(:,1)+s(:,m)-(canshu(2,m)-canshu(2,1))*C*b_joint_joint(:,1))];% h的前（M-1）行
            G1_joint_joint = [G1_joint_joint,2*(b_joint_joint(:,m)-b_joint_joint(:,1))];                                % G的前（M-1）列
            h1_jointDFT_joint = [h1_jointDFT_joint;(b_jointDFT_joint(:,m)-b_jointDFT_joint(:,1))'*(s(:,1)+s(:,m)-(canshu(6,m)-canshu(6,1))*C*b_jointDFT_joint(:,1))];% h的前（M-1）行
            G1_jointDFT_joint = [G1_jointDFT_joint,2*(b_jointDFT_joint(:,m)-b_jointDFT_joint(:,1))];  
            h1_sperate_joint = [h1_sperate_joint;(b_sperate_joint(:,m)-b_sperate_joint(:,1))'*(s(:,1)+s(:,m)-(canshu(4,m)-canshu(4,1))*C*b_sperate_joint(:,1))];% h的前（M-1）行
            G1_sperate_joint = [G1_sperate_joint,2*(b_sperate_joint(:,m)-b_sperate_joint(:,1))];                                % G的前（M-1）列
        end
        
        for m = 1:M
            h2_joint_joint = [h2_joint_joint;gm_joint_joint(:,m)'*s(:,m)];      % h的后M行
            G2_joint_joint = [G2_joint_joint,gm_joint_joint(:,m)];              % G的后M列
            h2_jointDFT_joint = [h2_jointDFT_joint;gm_jointDFT_joint(:,m)'*s(:,m)];      % h的后M行
            G2_jointDFT_joint = [G2_jointDFT_joint,gm_jointDFT_joint(:,m)];              % G的后M列
            h2_sperate_joint = [h2_sperate_joint;gm_sperate_joint(:,m)'*s(:,m)];      % h的后M行
            G2_sperate_joint = [G2_sperate_joint,gm_sperate_joint(:,m)];              % G的后M列
        end
        
        h_joint_joint = [h1_joint_joint;h2_joint_joint];
        G_joint_joint = [G1_joint_joint,G2_joint_joint];
        h_jointDFT_joint = [h1_jointDFT_joint;h2_jointDFT_joint];
        G_jointDFT_joint = [G1_jointDFT_joint,G2_jointDFT_joint];
        h_sperate_joint = [h1_sperate_joint;h2_sperate_joint];
        G_sperate_joint = [G1_sperate_joint,G2_sperate_joint];
        h_AOA = h2_sperate_joint;
        G_AOA = G2_sperate_joint;
        
        % 第一步普通最小二乘
        u_e_joint_joint = (G_joint_joint*G_joint_joint')\G_joint_joint*h_joint_joint;          % 第一步普通最小二乘定位结果
        err_joint_joint = norm(u_e_joint_joint - P); % 第一步普通最小二乘定位误差
        x_joint_joint=u_e_joint_joint(1);
        y_joint_joint=u_e_joint_joint(2);
        
        u_e_jointDFT_joint = (G_jointDFT_joint*G_jointDFT_joint')\G_jointDFT_joint*h_jointDFT_joint;          % 第一步普通最小二乘定位结果
        err_jointDFT_joint = norm(u_e_jointDFT_joint - P); % 第一步普通最小二乘定位误差
        x_jointDFT_joint=u_e_jointDFT_joint(1);
        y_jointDFT_joint=u_e_jointDFT_joint(2);
        
        u_e_sperate_joint = (G_sperate_joint*G_sperate_joint')\G_sperate_joint*h_sperate_joint;          % 第一步普通最小二乘定位结果
        err_sperate_joint = norm(u_e_sperate_joint - P); % 第一步普通最小二乘定位误差
        x_sperate_joint=u_e_sperate_joint(1);
        y_sperate_joint=u_e_sperate_joint(2);
        
        u_e_AOA = (G_AOA*G_AOA')\G_AOA*h_AOA;          % 第一步普通最小二乘定位结果
        err_AOA = norm(u_e_AOA - P); % 第一步普通最小二乘定位误差
        x_AOA=u_e_AOA(1);
        y_AOA=u_e_AOA(2);
%% 时延估计定位
        err_TOA=spotfi_Location_TOA(canshu(4,:),s,P);
%         if(TOA_flag~=[0 0 0])
%             err_sperate_joint=1000;
%             err_TOA=1000;
%         end
%% %%%%%%%%%%%%%%%%%%%%%%
Poserr=([err_joint_joint;err_sperate_joint;err_AOA;err_TOA;err_jointDFT_joint]).^2;
 end
 