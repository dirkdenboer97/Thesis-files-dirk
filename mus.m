function mu = mus(mu)

figure
plot(mu(1,:));
hold on
plot(mu(2,:));
legend('mu 1','mu 2');
ylabel('Value mu');
xlabel('Time steps')

end