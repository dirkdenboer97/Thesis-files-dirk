function [xf] = ukf_fcn_error(functions,z,u,Q,R,x0,P0,l,phi)


f = functions.statefcn;
h = functions.measfcn;

nx = length(x0);
nn = size(Q,2);
ny = size(R,1);

%Initialize vectors
xc = x0;
L  = nx ;      % Dimension of augmented state
ns = 2 * L + 1;         % Number of sigma points
Pf = P0;
R_s = R;

%Define UKF parameters
alpha = 0.1; 
beta = 2; 
kappa = 3; 
lambda = alpha ^2 * (L + kappa) - L;
gamma  = sqrt(L + lambda);
W_m_0  = lambda / (L + lambda);
W_c_0  = W_m_0 + 1 - alpha^2 + beta;
W_i    = 1/(2*(L + lambda));


%% Simulate true system
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
    
    %Create augmented system
    P_a     = Pf;
    x_a     = xc;
    
    %Generate sigma points
    [usvd, ssvd] = svd(P_a);
    gamma_sqrt_P_a = gamma * usvd * sqrt(ssvd);
    s_i = [x_a, ...
               repmat(x_a, 1, L) + gamma_sqrt_P_a, ...
               repmat(x_a, 1, L) - gamma_sqrt_P_a];
        
    %Prediction step
    for i=1:ns
    uk = u(:,k);
    X(:,i) = f(s_i(1:nx,i),u(:,k),zeros(15,1),Cbn(:,:,k));
    Z(:,i) = h(s_i(1:nx,i),zeros(ny,1));
    end
    
    % Expected prediction and measurement
    x_km = W_m_0 * X(:, 1) + W_i * sum(X(:, 2:end), 2);
    z_km = W_m_0 * Z(:, 1)   + W_i * sum(Z(:, 2:end), 2);

    % Remove expectations from X_x_kkm1 and Z_kkm1.
    X = bsxfun(@minus, X, x_km);
    Z   = bsxfun(@minus, Z, z_km);
    
    Inn_pp(:,k) = (znew(:,k)-z_km);
    Inn_pp2(:,:,k) = Inn_pp(:,k)*Inn_pp(:,k)';
    
    % Calculate covariance of the prediction.
    P_p =   (W_c_0 * X(:, 1)) * X(:, 1).' ...
             + W_i * (X(:, 2:end) * X(:, 2:end).') + QQ; 
    
    % Covariance of predicted observation
    P_zz =   (W_c_0 * Z(:, 1)) * Z(:, 1).' ...
           + W_i * (Z(:, 2:end) * Z(:, 2:end).') + R_s;
    
    % Covariance of predicted observation and predicted state
    P_xz =   (W_c_0 * X(:, 1)) * Z(:, 1).' ...
           + W_i * (X(:, 2:end) * Z(:, 2:end).');
       
    % Kalman gain
    K = P_xz / P_zz;
    
    % Correct the state.
    xf(:,k) = x_km + K * (z(1:ny,k) - z_km);
    xc = xf(:,k);
    
    %Residuals
    Inn_rr(:,k) = (znew(:,k)-h(xc,zeros(ny,1)));
    Inn_rr2(:,:,k) = Inn_rr(:,k)*Inn_rr(:,k)';
    
    % Correct the covariance.
    Pf = P_p - K * P_zz * K.';
    
    %IBAE 
    win = 100;
    Spp = zeros(4);
    if k > win
        Spp = sum(Inn_pp2(:,:,k-win+1:k),3)./win;
        H = zeros(4,15);H(1,1)=1;H(2,2)=1;H(3,4)=1;H(4,5)=1;
        R_s = Spp - H*(P_p*H');
        R_s = cov(Inn_pp(:,k-win+1:k)');
    end



end





end