function res = cbn(phi)

    phi_x = 0;
    phi_y = 0;
    phi_z = phi;
    
    res = [ cos(phi_y)*cos(phi_z)-sin(phi_x)*sin(phi_y)*sin(phi_z) ...
                    min(max(cos(phi_y)+sin(phi_x)*sin(phi_y)*sin(phi_z),-0.99),0.99)*sin(phi_z)...
                    -cos(phi_x)*sin(phi_y);
            min(max(-cos(phi_x)*sin(phi_z),-0.99),0.99) cos(phi_x)*cos(phi_z) sin(phi_x) ;
            sin(phi_y)*cos(phi_z)-sin(phi_x)*cos(phi_y)*sin(phi_z) ...
                    -sin(phi_x)*cos(phi_y)*cos(phi_z)+sin(phi_y)*sin(phi_z) cos(phi_x)*cos(phi_y)];
end