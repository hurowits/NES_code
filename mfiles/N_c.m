colors = 'rgbckym'
sigma = 1;
s = logspace(-4,2,1e5);

figure;
axes('FontSize',24);
grid on
hold on;

c=1;
for sigma = logspace(-1,1,7)
lambda=(log(sinh(sigma*2)./(sigma*2))-2*s);
plot(s(lambda>0), (1./(lambda(lambda>0))),colors(c),'LineWidth',4);
s_c=log(sinh(sigma*2)./sigma/2)/2;
plot([s_c s_c],[1e-2 1e4],'--k','LineWidth',2)
c=c+1;
end
xlabel('SMF');
ylabel('1/\lambda_c')
axtype(3);
axis(10.^([-3 2 -2 4]))
% legend([repmat('\sigma = ',7,1),num2str( logspace(-1,1,7)','%10.1e\n')])

%%
figure;
axes('FontSize',24);
grid on
hold on;

c=1;
for sigma = logspace(-1,1,7)
lambda=(log(sinh(sigma)./(sigma))-s);
plot(s(lambda>0), (1./(lambda(lambda>0))),colors(c),'LineWidth',4);
s_c=log(sinh(sigma)./sigma);
plot([s_c s_c],[1e-2 1e5],'--k','LineWidth',2)
c=c+1;
end
xlabel('SMF');
ylabel('1/\lambda')
axtype(3);
axis(10.^([-4 2 -2 5]))
% legend([repmat('\sigma = ',7,1),num2str( logspace(-1,1,7)','%10.1e\n')])
fileName = [DataPath,'largeN.eps'];
% print(gcf, '-depsc2', fileName);