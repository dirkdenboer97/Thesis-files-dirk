clear all
clc

%Localization 1 
dt = 0.01; %Sampling time
t = dt:dt:240; %4 minutes simulation time
l = length(t); %Number of time steps
u = zeros(2,l); %Initialze input vector
hc = 40; % 0.4 seconds of transition between inputs

%Generate inputs
sav = [min(0.004,abs(covdraw(1e-5))) -min(0.004,abs(covdraw(1e-4))) min(0.004,abs(covdraw(1e-5)))];
acv = [min(0.1,abs(covdraw(0.2))) min(0.05,abs(covdraw(0.1))) -min(0.1,abs(covdraw(0.2)))];
u(:,1:l/3) = [sav(1)*ones(l/3,1) ones(l/3,1).*acv(1)]';
u(:,l/3+1:l/3+hc) = [linspace(sav(1),sav(2),hc)' ones(hc,1).*acv(1)]';
u(:,l/3+hc+1:2*l/3) = [sav(2)*ones(l/3-hc,1) ones(l/3-hc,1).*acv(2)]';
u(:,2*l/3+1:2*l/3+hc) = [linspace(sav(2),sav(3),hc)' ones(hc,1).*acv(2)]';
u(:,2*l/3+1+hc:l) = [sav(3)*ones(l/3-hc,1) ones(l/3-hc,1).*acv(3)]'; 

%Define process- and measurement model
functions = struct('statefcn',@f6,'measfcn',@h6);

%Settings for  simulations
x0 = [0;0;0;15+(rand-0.5)*8;0;0;0;0;0;0]; %Variable initial condition long. velocity
R = eye(4).*[16 16 1 1]; %Measurement noise covariance matrix GNSS
Np = 50; %Number of particles
Npa = Np; 
Nt = Np/3; %Resampling threshold
Nta = Nt;

%Obtain measurements
[xa,xic,zins,zgps,zgps10hz,V,gmQ,Var] = actualsys_noise(x0,l,R,u);

%Extract measurement data
et = zins(1:2,:)-xa(1:2,:);
phi = xic(3,:);
l = l/10;
xa10hz = xa(:,1:10:end);
zins10hz = zins(:,1:10:end);
et10hz = et(:,1:10:end);
phi10hz = phi(:,1:10:end);
u = zins10hz;
z = [zins10hz(1:2,:);zins10hz(4:5,:)] - zgps10hz;

%Define uncertainties Q matrix
eta_acc = [1.17e-7;1.17e-7;1.17e-7];
eta_gyro = [1e-8;1e-8;1e-8];
eta_bacc = [2.1e-10;2.1e-10;2.1e-10];
eta_bgyro = [2.16e-8;2.16e-8;2.16e-8];
Q = diag([1e-30;1e-30;1e-30 ; eta_acc; eta_gyro; eta_bacc; eta_bgyro]);

%Define inital conditions filters
x0 = zeros(15,1);
P0 = abs(diag(covdraw(eye(15))));

%Simulate filters
functions = struct('statefcn',@errormodel2,'measfcn',@errormodelm);
functionsekf = struct('statefcn',@ekffunerror,'measfcn',@ekfmfunerror );
tic;[x_ckf,Pf_ckf] = ckf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi);disp('Time_ckf:');toc;
tic;[x_ukf] = ukf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi);disp('Time_ukf:');toc;
tic;[x_packf2,prop,Ppp,QQ,stanR,xp,znew,mu] = packf2_fcn_error_report(functions,z,u,Q,R,x0,P0,l,Np);disp('Time_packf2:');toc;
tic;[x_pfckf] = pfckf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi,Np,Nt);disp('Time_pfckf:');toc;
tic;[x_ckfapf] = ckfapf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi,Np,Nt);disp('Time_pfckf:');toc;

%Plots and tables
%mus(mu);
[ckf_err,ukf_err,packf2_err,pfckf_err,ckfapf_err] = rmsecalc(x_ckf,x_ukf,x_packf2,x_pfckf,x_ckfapf,zins10hz,l,xa10hz);
%plotfig(t,zgps10hz,xa,stanR,Ppp,xp,zins10hz)
rmsetab = rmse(ckf_err,ukf_err,packf2_err,pfckf_err,ckfapf_err);
filters = cell2table({'CKF','IBAE-UKF','PACKF','PF','CKFAPF'}');
rmsetab = [filters,array2table(rmsetab, 'VariableNames',{'Sector 1','Sector 2','Sector 3',...
    'Sector 4','Sector 5','Sector 6','Sector 7','Average',})]
%errorplot(ckf_err,ukf_err,pfckf_err,ckfapf_err,packf2_err,length(t)/10);
tracks(zgps10hz,zins10hz);
errorplot2(ckf_err,ukf_err,pfckf_err,ckfapf_err,packf2_err,length(t)/10);











