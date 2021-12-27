function [xpf] = ckfapf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi,Np,Nt)


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
win = 100;
inno = zeros(ny,win+l);
R_s = R;



%% Simulate system


for k=1:l

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


% Use newest observation    
    znew(:,k) = z(:,k);
    
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
    
    %% Update step CKF
    Xn=repmat(xp,1,m)+Sp(:,1:nx)*Cpt;
    for i = 1:2*nx
        Zp(:,i) = h(Xn(:,i),zeros(ny,1));
    end
    zhat=sum(Zp,2)/m;  
    inno(1:ny,k) = znew(1:ny,k)-zhat; 
    Zinn = 1/(sqrt(m))*(Zp-repmat(zhat,1,m));
    Rsr = chol(R);                                
    [~, r] = qr([Zinn Rsr]');
    Szz = r';
    Chinn = 1/(sqrt(m))*(Xn-repmat(xp,1,m));
    Pxz = Chinn*Zinn';
    K = (Pxz/Szz(:,1:nx)')/Szz(:,1:nx);
    obs = znew(:,k);
    xf(:,k) = xp + K*(obs-zhat);
    xd = xf(:,k);
    [~,r] = qr([Chinn-K*Zinn K*Rsr]'); 
    Sf = r';
    Pff(:,:,k) = Sf*Sf';
    pyoyo = Pff(:,:,k);
    
% Particle filter
    
    for i = 1:Np
      parn(:,i) = xd +  covdraw(pyoyo)   ; 
      ypar(:,i) = h(parn(:,i),zeros(ny,1));
      w(i,1) = (wp(i,k)^0)*max(1e-50,mvnpdf(znew(1:ny,k)-ypar(1:ny,i),zeros(ny,1), R_s(1:ny,1:ny)));
    end
    
    w = w./(sum(w)); 
    xpfo(:,k) = (w'*parn')';
    ypf(:,k) = h(xpfo(:,k),zeros(ny,1)); 
    
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
    else
        xpf(:,k) = xpfo(:,k);
    end
    xc = xpf(:,k);

    

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