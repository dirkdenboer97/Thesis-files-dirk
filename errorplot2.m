function errorplot2(ckf_err,ukf_err,pfckf_err,ckfapf_err,packf2_err,l)
figure
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
legend('Zero-line','CKF','PF','CKFAPF','IBAE-UKF','PACKF' ,'Interpreter','latex');%,'IBAE-UKF'
xlabel('Time steps')
ylabel('X-pos. error [m]')

figure
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
legend('Zero-line','CKF','PF','CKFAPF','IBAE-UKF','PACKF' ,'Interpreter','latex');%,'IBAE-UKF'
xlabel('Time steps')
ylabel('Y-pos. error [m]')




end


