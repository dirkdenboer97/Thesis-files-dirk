function [xnew,B] = errormodel2(x,u,noise,Cbn)

%Define INS measurements
L = u(1);
lambda = u(2);
h = u(3);
V_E = u(4);
V_N = u(5);
V_U = u(6);
f_E = u(7);
f_N = u(8);
f_U = u(9);

%Define constants
dt = 0.1; 
a_a = 0.005; 
a_g = 0.005; 
Re = 6370000; 
ff = 1/298.257;
Rn = Re*(1+ff*sin(L)^2);
Rm = Re*(1-2*ff+3*sin(L)^2);
g = 9.81; 
w_ie = 7.29e-5;

%Define matrices
Frr = [ 0 0 -V_N/(Rm+h)^2; V_E*sin(L)/((Rn+h)*cos(L)^2) 0 (-V_E)/((Rn+h)*cos(L)); 0 0 0];

Fvr = [ (2*w_ie*(V_U*sin(L)+V_N*cos(L)))+(V_E*V_U)/((Rn+h)*cos(L)^2),0,(V_E*V_U-V_E*V_N*tan(L))/(Rm+h)^2 ;
         -2*w_ie*V_E*cos(L)-(V_E^2)/((Rn+h)*cos(L)^2),0,(V_N*V_U-(V_E^2)*tan(L))/(Rm+h)^2 ;
         -2*w_ie*V_E*sin(L),0,(-V_E^2)/(Rn+h)-(V_N^2)/(Rm+h)^2+(2*g)/(Re+h) ];


Fvv = [ V_N*tan(L)/(Rm+h)-(V_U)/(Rn+h),2*w_ie*sin(L)+(V_E*tan(L))/(Rn+h),-2*w_ie*cos(L)-(V_E)/(Rn+h) ;
        -2*w_ie*sin(L)-(2*V_E*tan(L)/(Rn+h)),-V_U/(Rm+h),(-V_N/(Rm+h)) ;
        2*w_ie*cos(L)+(2*V_E)/(Rn+h),V_N/(Rm+h),0];
    
    
Fve = [ 0 f_U -f_N ; -f_U 0 f_E ; f_N -f_E 0];

Fer = [ 0,0,(-V_N)/(Rm+h)^2 ; w_ie*sin(L),0,(V_E)/(Rn+h)^2 ; 
        -w_ie*cos(L)-(V_E)/((Rn+h)*(cos(L)^2)),0,(V_E*tan(L))/(Rn+h)^2 ];
    
Fev = [ 0,1/(Rm+h),0 ; -1/(Rn+h),0,0 ; (-tan(L))/(Rm+h),0,0 ];


Fee = [0,w_ie*sin(L)+(V_E*tan(L))/Rn,-w_ie*cos(L)-(V_E/(Rn+h)) ;
        -w_ie*sin(L)-(V_E*tan(L)/Rn),0,-V_N/(Rm+h)  ;
        w_ie*cos(L)+V_E/(Rn+h),V_N/(Rm+h),0  ];


Fpp = [ 0,V_E*tan(L)/(Rm+h),-V_E/(Rn+h) ; 0,0,-V_N/(Rm+h) ; 0,0,0 ];

Fpv = eye(3);


Fvp = [0,2*w_ie*(V_U*sin(L)+V_N*cos(L))/(Rm+h)+(V_E*V_U)/((Rm+h)*(Rn+h)*cos(L)^2),V_E*(V_U-V_N*tan(L))/(Rn+h)^2;
       0,(-2*w_ie*V_E*cos(L)/(Rm+h))-(V_E^2)/((Rm+h)*(Rn+h)*cos(L)^2),...
       (-2*w_ie*V_E*cos(L))/(Rm+h)-(V_E^2)/((Rm+h)*(Rn+h)*cos(L)^2) ;
       0,(-2*w_ie*V_E*sin(L))/(Rm+h),(-V_N^2)/(Rm+h)^2-(V_E^2)/(Rn+h)^2+(2*g)/(Re+h) ];
   
Fep = [ 0,0,(-V_N)/(Rm+h)^2  ; 
        0,w_ie*sin(L)/(Rm+h),V_E/(Rn+h)^2 ; 
        0,(-w_ie*cos(L))/(Rm+h)-V_E/((Rm+h)*(Rn+h)*cos(L)^2),(V_E*tan(L))/(Rn+h)^2 ];


I = eye(3);
O = zeros(3);

A = [   I+Fpp*dt    Fpv*dt      O           O           O ;
        Fvp*dt      I+Fvv*dt    Fve*dt      -Cbn*dt      O ;
        Fep*dt      Fev*dt      I+Fee*dt    O       	-Cbn*dt;
        O           O           O           I-a_a*dt    O;
        O           O           O           O           I-a_g*dt];

B = [ O O O O ; Cbn O O O ; O Cbn O O ; O O I O ; O O O I];

xnew = A*x + noise;










        

end