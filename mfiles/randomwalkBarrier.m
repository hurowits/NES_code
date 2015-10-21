% numRealizations=1e4;
% N=1000;

%%
a = pi^2*N/2;
figure;
axes('FontSize',24);
hold on
plot((-a./B.^2),log((1:numRealizations)),'.','LineWidth',4)
grid on;
plot([-12 0],[2 14],'k','LineWidth',4)
xlabel('-\alpha/B^2')
ylabel('ln [Prob(Barrier<B)]')
% print(gcf, '-depsc2', 'CDF_B_sim.eps');

%%
sigmaU=sqrt(N);
x=linspace(0,200,1000);
h=hist(B,x);
figure;
axes('FontSize',24);

hold on
grid on
plot(x,h/sum(h)/(x(2)-x(1)),'r','LineWidth',4)
sigmaU=sqrt(N);
a=pi^2*sigmaU^2/2;

f=zeros(1,length(x));
for n=1:2:101
f = f+ 2*a./x.^3.*exp(-a*n^2./x.^2).*(2*a*n^2./x.^2-1)*8/pi^2;

end
% plot(x,2*a./x.^3.*exp(-a./x.^2).*(2*a./x.^2-1),'g','LineWidth',4);
plot(x,f,'--k','LineWidth',4);
plot(x,2*a./x.^3.*exp(-a./x.^2),'b','LineWidth',4)

legend(['simulation   ';...
        'expression   ';
        'approximation']);
    xlabel('R');
    ylabel('f(R)');
%     print(gcf, '-depsc2', 'PB_sim.eps');

%%

x1=linspace(0,1,1000);
x=linspace(0,1,1000);


f1 = log(1+x) - log(1-x);
f2 = (1-x.^2).^N;
f3 = ((1+x)./(1-x)).^(N*x);
p= f1./(f2.*f3);
% p2 = p(~isnan(p));
% p2=p2/sum(p2);

dx1=x1(2)-x1(1);
h=hist(maxR/N,x1);
figure;
hold on
plot(x,2*x.*exp(-N*x.^2)*N,'--k','LineWidth',4);
plot(x,p*N,'--r');
plot(x1,h/sum(h)/dx1);
axis([0 0.1 0 30])