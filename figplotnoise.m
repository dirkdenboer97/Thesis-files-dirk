function gis  = figplotnoise(z,znew,inno_win,V,k,win,innn)


subplot(1,2,1)
plot(V(1,1:10:end))
hold on
plot(innn(1,:))
legend('$\mathbf{v}$','$\mathbf{\tilde{y}}$','Interpreter','latex')
xlabel('Time steps')
ylabel('x-pos. deviation [m]')

subplot(1,2,2)
plot(V(1,1:10:end))
hold on
plot(innn(1,:))
legend('$\mathbf{v}$','$\mathbf{\tilde{y}}$','Interpreter','latex')
xlabel('Time steps')
ylabel('x-pos. deviation [m]')

% subplot(2,2,3)
% plot(z(1,:))
% hold on
% plot(znew(1,:))
% legend('$\mathbf{z}$','$\mathbf{z}_{pp}$','Interpreter','latex')
% xlabel('Time steps')
% ylabel('x-pos. [m]')
% 
% subplot(2,2,4)
% plot(z(1,:))
% hold on
% plot(znew(1,:))
% legend('$\mathbf{z}$','$\mathbf{z}_{pp}$','Interpreter','latex')
% xlabel('Time steps')
% ylabel('x-pos. [m]')



aa=2;





aa=6;
 m1 = moment(V(:,k-win:k),1,2)
 m2 = moment(inno_win,1,2)
 m3 = moment(V(:,k-win:k),2,2)
 m4 = moment(inno_win,2,2)
 m5 = moment(V(:,k-win:k),3,2)
 m6 = moment(inno_win,3,2)









end