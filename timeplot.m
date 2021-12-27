clear all
clc

figure
subplot(1,2,1)
X = categorical({'Standard PF','CKF pred.','KS-test','Our PF','CKF upd.','KDE','Resampling','Adaptive Q-matrix'});
X = reordercats(X,{'Standard PF','CKF pred.','KS-test','Our PF','CKF upd.','KDE','Resampling','Adaptive Q-matrix'});
Y = [0.732  0.676 0.654  0.392 0.266 0.247  0.057  0.052];

b = bar(X,Y,'FaceColor',[0.4660 0.6740 0.1880]);
xtips2 = b.XEndPoints;
ytips2 = b.YEndPoints;
labels2 = string(b.YData);
ylabel('Time [ms]')
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b.FaceColor = 'flat';
b.CData(1,:) = [0.9290 0.6940 0.1250];
b.CData(7,:) = [0.9290 0.6940 0.1250];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(5,:) = [0.8500 0.3250 0.0980];
%b.CData(5,:) = [0.4660 0.6740 0.1880];

subplot(1,2,2)
X = categorical({'PACKF','Standard PF','CKF'});
X = reordercats(X,{'PACKF','Standard PF','CKF'});
Y = [2.287 1.465 0.923];

b = bar(X,Y,'FaceColor',[0.3010 0.7450 0.9330]);
xtips2 = b.XEndPoints;
ytips2 = b.YEndPoints;
labels2 = string(b.YData);
ylabel('Time [ms]')
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b.FaceColor = 'flat';
b.CData(1,:) = [0.4660 0.6740 0.1880];
b.CData(2,:) = [0.9290 0.6940 0.1250];
b.CData(3,:) = [0.8500 0.3250 0.0980];

