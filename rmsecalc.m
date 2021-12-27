function [ckf_err,ukf_err,packf2_err,pfckf_err,ckfapf_err] = rmsecalc(x_ckf,x_ukf,x_packf2,x_pfckf,x_ckfapf,zins10hz,l,xa10hz)





ny=4;
for k = 1:l
    ckf_pos(:,k) = (zins10hz(1:ny,k)-x_ckf(1:ny,k));
    ukf_pos(:,k) = (zins10hz(1:ny,k)-x_ukf(1:ny,k));
    packf2_pos(:,k) = (zins10hz(1:ny,k)-x_packf2(1:ny,k));
    pfckf_pos(:,k) = (zins10hz(1:ny,k)-x_pfckf(1:ny,k));
    ckfapf_pos(:,k) = (zins10hz(1:ny,k)-x_ckfapf(1:ny,k));
    ckf_err(:,k) = xa10hz(1:ny,k)-ckf_pos(1:ny,k);
    ukf_err(:,k) = xa10hz(1:ny,k)-ukf_pos(1:ny,k);
    packf2_err(:,k) = xa10hz(1:ny,k)-packf2_pos(1:ny,k);
    pfckf_err(:,k) = xa10hz(1:ny,k)-pfckf_pos(1:ny,k);
    ckfapf_err(:,k) = xa10hz(1:ny,k)-ckfapf_pos(1:ny,k);
end


end