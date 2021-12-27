function trackk = tracks(zgps10hz,zins10hz)



figure
plot(zgps10hz(1,:),zgps10hz(2,:),zins10hz(1,:),zins10hz(2,:),...
    'LineWidth',2);%,z(1,:),z(2,:),z_imu(1,:),z_imu(2,:)) %zins10hz(1,:)-x_packf2(1,:),zins10hz(2,:)-x_packf2(2,:)
legend('GPS','INS','Interpreter','latex');
xlabel('x-position [m]');
ylabel('y-position [m]');


% figure
% subplot(1,3,1)
% plot(zgps10hz(1,:,1),zgps10hz(2,:,1),zins10hz(1,:,1),zins10hz(2,:,1),...
%     'LineWidth',2);%,z(1,:),z(2,:),z_imu(1,:),z_imu(2,:)) %zins10hz(1,:)-x_packf2(1,:),zins10hz(2,:)-x_packf2(2,:)
% legend('GPS','INS','Interpreter','latex');
% xlabel('x-position [m]');
% ylabel('y-position [m]');
% subplot(1,3,2)
% plot(zgps10hz(1,:,2),zgps10hz(2,:,2),zins10hz(1,:,2),zins10hz(2,:,2),...
%     'LineWidth',2);
% subplot(1,3,3)
% plot(zgps10hz(1,:,3),zgps10hz(2,:,3),zins10hz(1,:,3),zins10hz(2,:,3),...
%     'LineWidth',2);


end