function [x,f]=ploTcdf(rmse)
N=size(rmse,2);
f=zeros(size(rmse,1)+1,N);
x=zeros(size(rmse,1)+1,N);
for n=1:N
   [f(:,n),x(:,n)]=ecdf(rmse(:,n));       
end


end
