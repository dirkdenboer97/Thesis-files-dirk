function [xa,ximu,zins,zgps,zgps10hz,V,gmQ,Var] = actualsys_noise(x0,l,R,u)
xic = x0;
xta = x0;
xbc = x0;
nx = size(x0,1);
eps = [0;0;0;0];
Qn = eye(10).*[0;0;0;0;0;1e-8;1.17e-7;1.17e-7;2.1e-10;2.16e-8];

%Gaussian mixture
muQ = [-0.8 1.5; -1 0.4];
sigmaQ = cat(3,[0.9*sqrt(R(1,1)) 0.9*sqrt(R(2,2))],[0.1*(300.*sqrt(R(1,1))) 0.1*(300.*sqrt(R(2,2)))]); 
gmQ = gmdistribution(muQ,sigmaQ); 
muQ2 = [-0.1 0; -0.1 0.2];
sigmaQ2 = cat(3,[0.9*sqrt(R(3,3)) 0.9*sqrt(R(4,4))],[0.1*(300.*sqrt(R(3,3))) 0.1*(300.*sqrt(R(4,4)))]); 
gmQ2 = gmdistribution(muQ2,sigmaQ2); 
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gmQ,[x0 y0]),x,y);



for k=1:l
      
% Gaussian
     V(:,k) = covdraw(R);

% Gaussian mixture
if k >= 4000 && k<= 8000
    V(1:2,k) = random(gmQ);
    V(3:4,k) = random(gmQ2);
end

% Gaussian - varying variance
if k >= 8000 && k<= 12000
    V(:,k) = covdraw(max(0.2,rand)*5.*R);
end

% GPS outage
if k >= 16000 && k<= 17000
  V(:,k) = covdraw(100.*R);
end

% Random walk model
if k >= 17000
     V(:,k) = 0.6*eps + covdraw(R) ;
     eps = V(:,k);
end

% Gaussian
if k >= 21000
     V(:,k) = covdraw(R);
end

end

% Flicker noise
  V(1,12000:16000) = 100.*flicker(4001);
  V(2,12000:16000) = 100.*flicker(4001);
  V(3,12000:16000) = 10.*flicker(4001);
  V(4,12000:16000) = 10.*flicker(4001);

  
% Measurement data acquisition through noise injection
  for k=1:l
    uk = u(:,k);
    xa(:,k) = f6e(xta,zeros(nx,1),uk);
    xta = xa(:,k);  
    VV(:,k) = covdraw(Qn);
    [ximu(:,k),r_n(:,k),v_n(:,k),w_n(:,k),a_n(:,k),wb_n(:,k),ab_n(:,k)] = f6e(xbc,VV(:,k),uk);
    xbc = ximu(:,k);
    Cbn = cbn(xta(3));
    vve = Cbn*[xta(4:5);0];
    zgps(:,k) = [xta(1:2);vve(1:2)] + V(:,k);

end

zins = [r_n;v_n;a_n;w_n;ab_n;wb_n];
zgps10hz = zgps(:,1:10:end);
Var_postm = std(V(1,1:4000));
Var_postp = std(VV(1,1:4000));
Var = Var_postp*Var_postm;
mse = sqrt(mean((zgps10hz(1:2,:)-xa(1:2,1:10:end)).^2,2));
  

end











