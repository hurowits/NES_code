N=1000
% k=linspace(0,N,100);
%
% p = log((N+k)./(N-k))./(1-k.^2/N^2).^N;
%
% lnp = log(log((N+k)./(N-k))) - gammaln(N+k)-gammaln(N-k);

lnp = log(log((N+k)./(N-k))) - gammaln(N+k)-gammaln(N-k);

p1=exp(lnp/100);
C1=sum(p1);
p1=(p1/C1).^100
p1=p1./sum(p1);
figure;
plot(k/N,p1);

%%
% N=1000;
x=linspace(0,1,numRealizations/50);
dx = x(2)-x(1);
f1 = log(1+x) - log(1-x);
f2 = (1-x.^2).^N;
f3 = ((1+x)./(1-x)).^(N*x);
p= f1./(f2.*f3);
p2 = p(~isnan(p));
% p2=p2/sum(p2);
p3 = 2*x.*exp(-N*x.^2)*N;
% p3=p3./sum(p3);
sigma=sigma_vec;
x2=linspace(0,N/sigma,numRealizations/500);
dx=x2(2)-x2(1);


h=hist(maxV(:)/beta,x2*N/sigma); %sqrt(2)*beta is to normalize to unit steps
x = linspace(-1,160,200);
p3 = 2*x(x>0).*exp(-N*(x(x>0)).^2);

dx=x(2)-x(1);
h = hist(maxR,x);
figure;
axes('FontSize',24);
hold on;
grid on
% plot(x,p1,'--g','LineWidth',2);
% hist(maxV(:),x*N/sigma);
plot(x/N,h(1:end)/sum(h(1:end))/dx,'.r','MarkerSize',20);
% plot(x(~isnan(p)),p2,'b','LineWidth',6);
plot(x(x>0)/N,p3,'--g','LineWidth',4);
xlabel('u =U {\sigma}/{N}','Interpreter','tex')
legend(['hist(U_{max})';...
        'p(U_{max})   ']);
axis([0 1 0 0.1])
% print(gcf, '-depsc2', 'PmaxV');

%%
N=1000;
x=linspace(0,0.2,numRealizations/100);
fB = 2*exp(-N*x.^2).*(sqrt(pi*N/2) .* exp(x.^2*N/2).*(x.^2*N-1).*erf(x*sqrt(N/2))+N*x);
fB=fB./sum(fB);
h=hist(Barrier(:),x*sqrt(2*N)); %sqrt(2)*beta is to normalize to unit steps
figure;
axes('FontSize',24);
hold on;
grid on
% plot(x,p1,'--g','LineWidth',2);
% hist(maxV(:),x*N/sigma);
plot(x,h(1:end)/sum(h(1:end)),'.r','MarkerSize',20);
plot(x,fB,'b','LineWidth',6);
xlabel('b = B (2N)^{-1/2}','Interpreter','tex')
legend(['hist(b)      ';...
        'p(b) (Eq. 60)']);
% axis([0 0.1 0 0.018])
% print(gcf, '-depsc2', 'PB');

%%
X = [zeros(numRealizations,1) randn(numRealizations,N-1)];
R= cumsum(X,2);
%%
maxR =  (max(R,[],2));
minR =  (min(R,[],2));
B =   max(R,[],2)-min(R,[],2);

R2 = R(abs(R(:,end))<2,:);
maxR2 =  (max(R2,[],2));
minR2 =  (min(R2,[],2));
B2 =   max(R2,[],2)-min(R2,[],2);

[hMin, xMin] = hist(minR,100);
[hMax, xMax] = hist(maxR,100);

[hMin2, xMin2] = hist(minR2,100);
[hMax2, xMax2] = hist(maxR2,100);

figure;
axes('FontSize',24);

hold on

% plot(maxR,minR,'og');
plot(maxR2,minR2,'.r');
% 
% plot(xMax2,-200*hMax2/sum(hMax2),'--r','LineWidth',4)
% plot(200*hMin2/sum(hMin2),xMin2,'--b','LineWidth',4)
grid on
xlabel('x_{max}');
ylabel('x_{min}');
% legend(['(x_{max},x_{min})';...
%         'p(u^-)   ';...
%         'p(u^+)   '], 'Location','SouthEast')
% print(gcf, '-depsc2', 'MinMax');

figure;
axes('FontSize',24);

hold on


plot(maxR,minR,'.g');

plot(xMax,-100*hMax/sum(hMax),'--r','LineWidth',4)
plot(100*hMin/sum(hMin),xMin,'--b','LineWidth',4)
grid on
xlabel('u^{+}');
ylabel('u^{-}');
legend(['(u^+,u^-)';...
    'p(u^-)   ';...
    'p(u^+)   '], 'Location','SouthEast')
% axis([0 4 -4 0])
%%
x=linspace(1e-8,1,numRealizations/100);


p3 = 2*x.*exp(-N*x.^2)*N;
% p3=p3./sum(p3);
sigma=sigma_vec;
h=hist(maxV(:)/sqrt(2*N),x); %sqrt(2)*beta is to normalize to unit steps
h2=hist(sqrt(2)*maxR2/N,x);

figure;
axes('FontSize',24);
hold on;
grid on
% plot(x,p1,'--g','LineWidth',2);
% hist(maxV(:),x*N/sigma);
plot(x(1:end)*N,h(1:end)/sum(h(1:end))/((x(2)-x(1))),'ob','MarkerSize',10);
% plot(x(~isnan(p))*N,p2,'b','LineWidth',6);
plot(x*N,h2/sum(h2(1:end))/(x(2)-x(1)),'.r','MarkerSize',30);
plot(x*N,p3,'k','LineWidth',4);
% xlabel('\u =U (2N)^{-1/2}','Interpreter','tex')
xlabel('K');
ylabel('f(K)')
legend(['max[U(x)]  ';...
        'random walk';...
        'expression ']);
axis([0 0.1*N 0 31])
print(gcf, '-depsc2', 'f_K');