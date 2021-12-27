function [xnew,r_n,v_n,w_n,a_n,wb_n,ab_n] = f6e(x,w,u)
dt = 0.01;    
xs = x(1);
y = x(2);
phi = x(3);
vx = x(4);
vy = x(5);
wz = x(6);
ax = x(7);
ay = x(8);
wb = x(9);
ab = x(10);

delta = u(1);

Caf = 50000;
Car = 50000;
m = 1700;
lf = 1.5;
lr = 1.5;
Iz = ((0.5*(lf+lr))^2*m);

b1 = -2*(Caf+Car)/m;
b2 = 2*(-Caf*lf + Car*lr)/m;
b3 = 2*Caf/m;
b4 = 2*(-Caf*lf + Car*lr)/Iz;
b5 = -2*(Caf*lf^2 + Car*lr^2)/m;
b6 = 2*Caf*lf/Iz;

a_a = 0.0005; %TURNING KNOB default 1
a_g = 0.0005; %TURNING KNOB default 1
% taua = 1/a_a;
% taug = 1/a_g;
% 
% xnew = [    xs + vx*dt*cos(phi) - dt*vy*sin(phi) ;
%             y + vx*dt*sin(phi) + dt*vy*cos(phi) ;
%             phi + dt*(wz -  wb)            ;
%             vx+dt*(ax )                      ;
%             vy+dt*(ay  - ab)                      ;
%             wz + dt*((b4*vy/vx) + b5*wz/vx + b6*delta) ;
%             u(2) ;
%             (b1*vy/vx) + (b2/vx-vx)*wz + b3*delta;
%             wb/taug;
%             ab/taua  ]  + w ;
        
        xnew = [    xs + vx*dt*cos(phi) - dt*vy*sin(phi) ;
            y + vx*dt*sin(phi) + dt*vy*cos(phi) ;
            phi + dt*(wz -  wb)            ;
            vx+dt*(ax )                      ;
            vy+dt*(ay  - ab)                      ;
            wz + dt*((b4*vy/vx) + b5*wz/vx + b6*delta) ;
            u(2) ;
            (b1*vy/vx) + (b2/vx-vx)*wz + b3*delta;
            wb + dt*(-a_g*wb);
            ab + dt*(-a_a*ab)  ]  + w ;
        
phi_x = 0;
phi_y = 0;
phi_z = xnew(3);

Cbn = [ cos(phi_y)*cos(phi_z)-sin(phi_x)*sin(phi_y)*sin(phi_z) ...
                cos(phi_y)*sin(phi_z)+sin(phi_x)*sin(phi_y)*sin(phi_z) -cos(phi_x)*sin(phi_y);
        -cos(phi_x)*sin(phi_z) cos(phi_x)*cos(phi_z) sin(phi_x) ;
        sin(phi_y)*cos(phi_z)-sin(phi_x)*cos(phi_y)*sin(phi_z) ...
                -sin(phi_x)*cos(phi_y)*cos(phi_z)+sin(phi_y)*sin(phi_z) cos(phi_x)*cos(phi_y)];
           
r_n = [xnew(1);xnew(2);0];%MAYBE Cbn GONE
v_n = Cbn*[xnew(4);xnew(5);0];
w_n = Cbn*[0;0;xnew(6)];
a_n = Cbn*[xnew(7);xnew(8);0];
wb_n = Cbn*[0;0;xnew(9)];
ab_n = Cbn*[0;xnew(10);0];

end









