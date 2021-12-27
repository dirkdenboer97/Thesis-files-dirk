function [xf,prop,Ppp,QQ,stanR,xp,znew,mu] = packf2_fcn_error_report(functions,z,u,Q,R,x0,P0,l,Np)

%Extract process- and measurement functions
f = functions.statefcn;
h = functions.measfcn;

%Number of states and number of observations
nx = length(x0);
ny = size(R,1);

%Initialze matrices, vectors, variables
w = ones(Np,1).*(1/Np);
xf = zeros(nx,l);
xf(:,1) = x0;
xc = x0;
m = 2*nx;
Sf = chol(P0);
Cpt = sqrt(nx)*[eye(nx) -eye(nx)];
stanR = repmat(diag(R),1,l);
QQ = zeros(nx,nx,l);
Ppp = zeros(nx,nx,l);
prop = zeros(nx,Np,l);
win = 50;
binsize = 100;
mu = zeros(ny,l+1);
Gaussvec = zeros(ny,win+l);
R_s = R;
binss = ones(binsize,ny);
pd = ones(binsize,1);
Cbn = zeros(3,3,l);
Xp = zeros(nx,m);
Zp = zeros(ny,m);
znew = zeros(ny,l);
zhat = zeros(ny,l);
ypar = zeros(1,ny,Np);
ypfo = zeros(ny,l);
ypf = zeros(ny,l);
kap_pos = 50;
kap_vel = 5;
eta = 0.005;
rho = 0.01;
bw_par = 3.5; 
gamma = 0.1;


%% Simulate system

for k=1:l
%% Construct Q and R

 QQ(:,:,k) = eye(15).*diag(Q);
 if k > win+1
    minoo = mean(inno_win(1:ny,:),2);                                       %Mean of innovation sequence
    Qt = rho.*minoo.^2 + eta.*eye(4);                                       %Adjusted Q matrix
    QQ(1:2,1:2,k) = min(eye(2).*Qt(1:2),eye(2).*kap_pos);                   %Upper bound on position innovation
    QQ(4:5,4:5,k) = min(eye(2).*Qt(3:4),eye(2).*kap_vel);                   %Upper bound on velocity innovation
    R_s = R+1e10.*diag(max(sign(sum(minoo,2)...
        - [kap_pos;kap_pos;kap_vel;kap_vel]),0));                           %Increase R in case of GNSS outage
  end
    
   
%% Prediction step CKF
    
  X=repmat(xc,1,m)+Sf(:,1:nx)*Cpt;                                          %Draw cubature points
  for i = 1:2*nx
      [Xp(:,i),~] = f(X(:,i),u(:,k),zeros(15,1),Cbn(:,:,k));                %Propagate cubature points
  end
  xp=sum(Xp,2)/m;                                                         
  Chi = 1/(sqrt(m))*(Xp-repmat(xp,1,m));
  Qsr = chol(QQ(:,:,k));                                                  
  [~, r] = qr([Chi Qsr]');
  Sp = r';                                                                  %Square root of state covariance
  Ppp(:,:,k) = Sp*Sp';                                                      %State covariance
    
%% Particle filter
    
  znew(:,k) = z(:,k);

  %Draw particles
  prop(:,:,k) = covdraw(Ppp(:,:,k),Np);                                  
  parn = repmat(xp,1,Np)  +  prop(:,:,k);                                   %Sample particles around prediction
  ypar(:,:,:) = [parn(1:2,:);parn(4:5,:)];                                  %Extract measurable states

  %Weight particles
  if k > win+1
     m3d = repmat(repmat(znew(:,k)',1,1,Np)-ypar,binsize,1,1);              %Distance between observation and particle
     [~,index] = min(abs(repmat(binss(:,:),1,1,Np)-m3d),[],1);              %Index of measurement likelihood
     w(:,1) = max(1e-200,prod(pd(index(:,:,:)),2));                         %Calculate weight
  else
      for i = 1:Np
        w(i,1) = mvnpdf((znew(1:ny,k)'-ypar(:,:,i))'...                     %Standard particle weighting
            ,zeros(ny,1), R_s(1:ny,1:ny));
      end
  end

  wr = w./(sum(w));                                                         %Normalize particle weights
  ypfo(:,k) = h((wr'*parn')',zeros(ny,1));                                  %Particle filter observation

  %Apply Savitzky-Golay filter                                              %Filter over last 5 outputs            
  tv = 5;
  if k>tv
       smy = smoothdata(ypfo(:,k-tv:k),2,'sgolay',[tv 0]);
       ypf(:,k) = smy(:,end);
  else
       ypf(:,k) = ypfo(:,k);
  end

    
    %% Update step CKF

    Xn=repmat(xp,1,m)+Sp(:,1:nx)*Cpt;                                       %Draw new cubature points
    for i = 1:2*nx
        Zp(:,i) = h(Xn(:,i),zeros(ny,1));                                   %Propagate cubature points
    end
    zhat(:,k)=sum(Zp,2)/m;                                                  %Predicted observation
    Zinn = 1/(sqrt(m))*(Zp-repmat(zhat(:,k),1,m));
    Rsr = chol((ones(ny,1)-mu(:,k)).*R_s);                                  %Adjusted measurement covariance                             
    [~, r] = qr([Zinn Rsr]');
    Szz = r';
    Chinn = 1/(sqrt(m))*(Xn-repmat(xp,1,m));
    Pxz = Chinn*Zinn';
    K = (Pxz/Szz(:,1:nx)')/Szz(:,1:nx);                                     %Compute Kalman gain
    obs = (ones(ny,1)-mu(:,k)).*znew(:,k) + mu(:,k).*mean(ypf(:,k:k),2);    %Adjusted observation   
    xf(:,k) = xp + K*(obs-zhat(:,k));                                       %Final state estimate
    xc = xf(:,k);
    [~,r] = qr([Chinn-K*Zinn K*chol(R)]');                                  
    Sf = r';                                                                %Square root of state covariance


%% Check for Gaussianity and make distribution
    if k>win 
       inno_win(1:ny,:) = -(z(:,k-win:k)-zhat(:,k-win:k));                  %Innovations vector
       [binss,pd,stanR(:,k)] = makedis(inno_win,bw_par,win,binsize);        %Estimate measurement noise distribution
       [mu(:,k+1),Gaussvec(:,k+1)] = isgauss(inno_win,Gaussvec(:,k-win:k)); %Gaussianity check
       mu(:,k+1) = min(mu(:,k+1),1-gamma);                                  %Upper bound on Gaussianity parameter                           
    end

end



   


end