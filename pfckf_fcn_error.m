function [xpf] = pfckf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi,Np,Nt)


f = functions.statefcn;
h = functions.measfcn;

nx = length(x0);
ny = size(R,1);
nn = size(Q,2);
wp = ones(Np,1).*(1/Np);
w = ones(Np,1).*(1/Np);
parn = repmat(x0,1,Np);
xf = zeros(nx,l);
xf(:,1) = x0;
xc = x0;
m = 2*nx;
Sf = chol(P0);
Cpt = sqrt(nx)*[eye(nx) -eye(nx)];
win = 50;
mu = zeros(ny,l+1);
inno = zeros(ny,win+l);
R_s = R;



%% Simulate system


for k=1:l

% Use newest observation    
    znew(:,k) = z(:,k);
    
% Transformation matrix
    phi_x = 0;
    phi_y = 0;
    phi_z = phi(k);
    Cbn(:,:,k) = [ cos(phi_y)*cos(phi_z)-sin(phi_x)*sin(phi_y)*sin(phi_z) ...
                    cos(phi_y)*sin(phi_z)+sin(phi_x)*sin(phi_y)*sin(phi_z)...
                    -cos(phi_x)*sin(phi_y);
            -cos(phi_x)*sin(phi_z) cos(phi_x)*cos(phi_z) sin(phi_x) ;
            sin(phi_y)*cos(phi_z)-sin(phi_x)*cos(phi_y)*sin(phi_z) ...
                    -sin(phi_x)*cos(phi_y)*cos(phi_z)+sin(phi_y)*sin(phi_z) cos(phi_x)*cos(phi_y)] ;
% Q matrix    
    QQ = eye(15).*diag(Q);

    
% Prediction step CKF
    X=repmat(xc,1,m)+Sf(:,1:nx)*Cpt;
    for i = 1:2*nx
        [Xp(:,i),B] = f(X(:,i),u(:,k),zeros(15,1),Cbn(:,:,k)); 
    end
    xp=sum(Xp,2)/m;
    Chi = 1/(sqrt(m))*(Xp-repmat(xp,1,m));
    Qsr = chol(QQ);
    [~, r] = qr([Chi Qsr]');
    Sp = r';
    Ppp = Sp*Sp';
    
% Particle filter
    
    for i = 1:Np
      parn(:,i) = xp +  covdraw(Ppp)  ; 
      ypar(:,i) = h(parn(:,i),zeros(ny,1));
      w(i,1) = (wp(i,1)^0)*max(1e-200,mvnpdf(znew(1:ny,k)-ypar(1:ny,i),zeros(ny,1), R_s(1:ny,1:ny)));
    end

    w = w./(sum(w)); 
    xpfo(:,k) = (w'*parn')';

    % Resampling
    Neff = 1/sum(w.^2);
    if Neff < Nt
       [parn, w] = resample(parn, w, 'multinomial_resampling');
    end
    wp(:,k+1) = w;
    paro = parn;
    
    % Filter the PF output
    tv = 10;
    if k>tv
       smy = smoothdata(xpfo(:,k-tv:k),2,'sgolay',[tv 0]);
       xpf(:,k) = smy(:,end);
        aa = 5;
    else
        xpf(:,k) = xpfo(:,k);
    end
   
   xc = xpf(:,k);
   ypf(:,k) = h(xpf(:,k),zeros(ny,1));

% Update step CKF
    Xn=repmat(xp,1,m)+Sp(:,1:nx)*Cpt;
    for i = 1:2*nx
        Zp(:,i) = h(Xn(:,i),zeros(ny,1));
    end
    zhat=sum(Zp,2)/m;  
    inno(1:ny,k) = znew(1:ny,k)-zhat; 
    Zinn = 1/(sqrt(m))*(Zp-repmat(zhat,1,m));
    Rsr = chol((ones(ny,1)-mu(:,k)).*R);                        
    [~, r] = qr([Zinn Rsr]');
    Szz = r';
    Chinn = 1/(sqrt(m))*(Xn-repmat(xp,1,m));
    Pxz = Chinn*Zinn';
    K = (Pxz/Szz(:,1:nx)')/Szz(:,1:nx);
    obs = (ones(ny,1)-mu(:,k)).*znew(:,k) + mu(:,k).*mean(ypf(:,k:k),2);  
    xf(:,k) = xp + K*(obs-zhat);
    xc = xf(:,k);
    [~,r] = qr([Chinn-K*Zinn K*chol(R)]'); 
    Sf = r';
    Pff(:,:,k) = Sf*Sf';
    pyoyo = Pff(:,:,k);
    

% Completely use PF output
    mu(:,k+1) = ones(ny,1).*0.9999;
       
    

end




function [xk, wk, idx] = resample(xk, wk, resampling_strategy)

Ns = length(wk);  
switch resampling_strategy
   case 'multinomial_resampling'
      with_replacement = true;
      idx = randsample(1:Ns, Ns, with_replacement, wk);

   case 'systematic_resampling'
      edges = min([0 cumsum(wk)'],1); 
      edges(end) = 1;                
      u1 = rand/Ns;
      [~, idx] = histc(u1:1/Ns:1, edges);
   otherwise
      disp('Resampling strategy not implemented')
end

xk = xk(:,idx);                    
wk = repmat(1/Ns, 1, Ns);  
wk = wk';

end  





end