function [H] = kstest2(x)


alpha  =  0.1;
x = x(~isnan(x));
n = length(x);
[sampleCDF,x] = ecdf(x);
x = x(2:end);
yCDF = normcdf(x,0,1);
nullCDF  =  yCDF;       
delta1    =  sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
deltaCDF  =  abs([delta1 ; delta2]);
KSstatistic   =  max(deltaCDF);
t = n * KSstatistic;
k = ceil(t):n;       % k is the variable of summation
pValue = sum(exp( log(t) - n*log(n)+ gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
k.*log(k - t) + (n-k-1).*log(t+n-k)));
 
%[~,ind] = max(delta1);

% figure
% plot(x,yCDF)
% hold on
% plot(x,sampleCDF(2:end))
% xline(x(ind),'-',{'$D_k$'},'Interpreter','latex')
% plot(x(ind),yCDF(ind),'ro')
% plot(x(ind),sampleCDF(ind),'ro')
% xlabel('Samples')
% ylabel('Cumulative probability')


H  =  (pValue < alpha);

end