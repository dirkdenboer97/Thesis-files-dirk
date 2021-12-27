function errorplot(ckf_err,ukf_err,pfckf_err,ckfapf_err,packf2_err,l)

% 
% ckf_err(1,:) = smoothdata(ckf_err(1,:),'sgolay',50);
% ukf_err(1,:) = smoothdata(ukf_err(1,:),'sgolay',50);
% pfckf_err(1,:) = smoothdata(pfckf_err(1,:),'sgolay',50);
% ckfapf_err(1,:) = smoothdata(ckfapf_err(1,:),'sgolay',50);
% packf2_err(1,:) = smoothdata(packf2_err(1,:),'sgolay',50);
% ckf_err(2,:) = smoothdata(ckf_err(2,:),'sgolay',50);
% ukf_err(2,:) = smoothdata(ukf_err(2,:),'sgolay',50);
% pfckf_err(2,:) = smoothdata(pfckf_err(2,:),'sgolay',50);
% ckfapf_err(2,:) = smoothdata(ckfapf_err(2,:),'sgolay',50);
% packf2_err(2,:) = smoothdata(packf2_err(2,:),'sgolay',50);

%X-pos
figure
subplot(2,3,1)
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

subplot(2,3,2)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);

subplot(2,3,3)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);


subplot(2,3,4)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);
xlabel('Time steps')
ylabel('X-pos. error [m]')


subplot(2,3,5)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);

subplot(2,3,6)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(1,:),'LineWidth',2);
plot(pfckf_err(1,:),'LineWidth',2);
plot(ckfapf_err(1,:),'LineWidth',2);
plot(packf2_err(1,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(1,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);



%Y-pos
figure
subplot(2,3,1)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
legend('Zero-line','CKF','PF','CKFAPF','IBAE-UKF','PACKF' ,'Interpreter','latex');%,'IBAE-UKF'
xlabel('Time steps')
ylabel('Y-pos. error[m]')

subplot(2,3,2)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);

subplot(2,3,3)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);

subplot(2,3,4)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);
xlabel('Time steps')
ylabel('Y-pos. error [m]')


subplot(2,3,5)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);

subplot(2,3,6)
plot(zeros(1,l),'LineWidth',2);
hold on
plot(ckf_err(2,:),'LineWidth',2);
plot(pfckf_err(2,:),'LineWidth',2);
plot(ckfapf_err(2,:),'LineWidth',2);
plot(packf2_err(2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
plot(ukf_err(2,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2);











end







