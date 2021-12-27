function [xf,Pff] = ckf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi)


f = functions.statefcn;
h = functions.measfcn;

nx = length(x0);
nn = size(Q,2);
ny = size(R,1);

%Draw cubature points
xf = zeros(nx,l);
xf(:,1) = x0;
xc = x0;
m = 2*nx;
Sf = chol(P0);
Cpt = sqrt(nx)*[eye(nx) -eye(nx)];
R_s = R;

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

% Prediciton step
    X=repmat(xc,1,m)+Sf(:,1:nx)*Cpt;
    for i = 1:2*nx
        [Xp(:,i),B] = f(X(:,i),u(:,k),zeros(15,1),Cbn(:,:,k)); 
    end
    xp=sum(Xp,2)/m;
    Chi = 1/(sqrt(m))*(Xp-repmat(xp,1,m));
    Qsr = chol(QQ);
    [~, r] = qr([Chi Qsr]');
    Sp = r';
  
% Update step
    Xn=repmat(xp,1,m)+Sp(:,1:nx)*Cpt;
    for i = 1:2*nx
        Zp(:,i) = h(Xn(:,i),zeros(ny,1));
    end
    zhat=sum(Zp,2)/m;  
    Zinn = 1/(sqrt(m))*(Zp-repmat(zhat,1,m));
    Rsr = chol(R_s);
    [~, r] = qr([Zinn Rsr]');
    Szz = r';
    Chinn = 1/(sqrt(m))*(Xn-repmat(xp,1,m));
    Pxz = Chinn*Zinn';
    K = (Pxz/Szz(:,1:nx)')/Szz(:,1:nx);
    Inn_pp(:,k) = (znew(1:ny,k)-zhat);
    xf(:,k) = xp + K*Inn_pp(:,k);
    xc = xf(:,k);
    [~,r] = qr([Chinn-K*Zinn K*Rsr]');
    Sf = r';
    Pff = Sf*Sf';
   
    

end




end