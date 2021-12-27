function plotfig(t,zgps10hz,xa,stanR,Ppp,xp,zins10hz)


xconf = [t(1,1:10:end) t(1,end:-10:1)] ;         
yconf = [zgps10hz(1,1:1:end)-3.*stanR(1,1:1:end)...
    zgps10hz(1,end:-1:1)+3.*stanR(1,end:-1:1)];
Ppdata(1,:) = Ppp(1,1,:);
yconf2 = [(zins10hz(1,1:1:end)-xp(1,1:1:end))-3.*Ppdata(1,1:1:end)...
    (zins10hz(1,end:-1:1)-xp(1,end:-1:1))+3.*Ppdata(1,end:-1:1)];

figure
p = fill(100.*xconf,yconf,'red');
p.FaceColor = [0.9290 0.6940 0.1250]; 
p.FaceAlpha = 0.3;
hold on
p2 = fill(100.*xconf,yconf2,'blue');
p2.FaceColor = [0.8 0.8 1]; 
p2.FaceAlpha = 0.5;
plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
plot(xa(1,:));
legend('$3\sigma$ obs. likelihood', '$3\sigma$ prop. distr.', 'GPS x-pos.','Actual x-pos.','Interpreter', 'latex')
xlabel('Time steps')
ylabel('x-position')
% 
% figure
% subplot(2,3,1)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% %plot(zgps(1,1:10:end),'r-');
% %plot(xp(1,1:10:end),'r-');
% %plot(zgps(1,1:1:end),'r-');
% %plot(xa(1,:)+V(1,:),'r-');
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));
% legend('$3\sigma$ obs. likelihood', '$3\sigma$ prop. distr.', 'GPS x-pos.','Actual x-pos.','Interpreter', 'latex')
% xlabel('Time steps')
% ylabel('x-position')
% 
% subplot(2,3,2)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));
% 
% subplot(2,3,3)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));
% 
% subplot(2,3,4)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));
% 
% subplot(2,3,5)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));
% 
% subplot(2,3,6)
% p = fill(100.*xconf,yconf,'red');
% p.FaceColor = [0.9290 0.6940 0.1250]; 
% p.FaceAlpha = 0.3;
% hold on
% p2 = fill(100.*xconf,yconf2,'blue');
% p2.FaceColor = [0.8 0.8 1]; 
% p2.FaceAlpha = 0.5;
% plot(100.*t(1:10:end),zgps10hz(1,1:1:end),'LineWidth',2);
% plot(xa(1,:));


end