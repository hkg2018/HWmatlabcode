

for idx=1:length(Res_tau0)
    for ii=1:length(SNR)
        load(['..\mat\' 'performance_all_estimator' num2str(Res_aoa/pi*180) 'deg'  num2str(Res_tau0(idx)*1e9) 'ns_' 'SNR' num2str(SNR(ii)) 'dB.mat']);
        
        OtherRmse=zeros(size(Pos_rmse,1)*size(Pos_rmse,2),5);
        for jj=1:size(OtherRmse,2)
            rmse0=Pos_rmse(:,:,jj);
            OtherRmse(:,jj)=rmse0(:);
        end
        [x,f]=ploTcdf([MusiCrmse(:),FastHRrmse(:),OtherRmse,MFrmse(:)]);
        
        x=x(:,[1,3,2,8,4,6,5,7]);
        f=f(:,[1,3,2,8,4,6,5,7]);
        
        figure;
        plot(x(:,3:7),f(:,3:7)*100,'LineWidth',2);
        grid on;
        xlabel('Position error/m');
        ylabel('percent/%');
        axis([x(1,1),2,0,100]);
%         title({[ 'ResTau=' num2str(Res_tau0(idx)*1e9) 'ns,' 'ResAoA=' num2str(Res_aoa/pi*180) 'deg']});
        legend('improved dpd-MF','dpd-MF','Joint-solution','TOA-music','AOA-music','location','best');
%         set(gca,'xtick',[.1,.2,.3]);
    end
end