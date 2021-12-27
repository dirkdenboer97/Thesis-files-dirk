function [xnew] = f6(x,w,u,z)

dt = 0.01;    
xs = x(1);
y = x(2);
phi = x(3);
vx = x(4);
vy = x(5);
wz = x(6);
ax = u(2);
ay = x(8);
wb = x(9);
ab = x(10);




xnew = [    xs + vx*dt*cos(phi) - dt*vy*sin(phi) ;
            y + vx*dt*sin(phi) + dt*vy*cos(phi) ;
            phi + dt*(wz -  wb)            ;
            vx+dt*(ax )                      ;
            vy+dt*(ay  - ab)                      ;
            z(3) ;
            z(4);
            z(5);
            wb
            ab  ]  + w ;
        
        
        
     
end


