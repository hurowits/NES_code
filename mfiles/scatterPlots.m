% CDFFigures(CurrentEnsemble,'Current',1)
% CDFFigures(abs(CurrentEnsemble),'|Current|',1)
% CDFFigures(log10(abs(CurrentEnsemble)),'log10|Current|',2)
%% Scatter plots
% figure;
% axes('FontSize',24);
% hold on;
% grid on;
% plot(barrier,abs(CurrentEnsemble),'.','MarkerSize',0.1)
% xlabel('Barrier height');
% ylabel('|Current|');
% 
% figure;
% axes('FontSize',24);
% hold on;
% grid on;
% plot(barrier,log10(abs(CurrentEstimateEnsemble)),'.','MarkerSize',0.1)
% xlabel('Barrier height');
% ylabel('log10(|Current|)');
% 
% figure;
% axes('FontSize',24);
% hold on;
% grid on;
% plot(barrier,log10(abs(CurrentEnsemble)),'.','MarkerSize',0.1)
% xlabel('Barrier height');
% ylabel('log10(|CurrentEstimate|)');
%%
flagPrint=2;
data = barrier;

xB = linspace(min(barrier),max(barrier),100);
figure;
axes('FontSize',24);
hold on;
grid on
% xB=1:150;
p = polyfit(sort(data)',log((length(data):-1:1)/length(data)),2);
muB = -p(2)/(2*p(1));
sigmaB = sqrt(1/2/abs(p(1)))

%         p3 = polyfit(sort(data)',log10((1:length(data))/length(data)),2);
% plot(xB,log(1-0.5*(erf((xB-muB)/sqrt(2*sigmaB^2))+erf((xB+muB)/sqrt(2*sigmaB^2))))+log(2*sqrt(pi))-2*log(p(3)),'--b','LineWidth',5)
% plot(xB,log(1-erf(xB/sqrt(2*sigmaB^2)))+2*log(2*sqrt(pi))+p(3),'--b','LineWidth',5)

% plot(x,log(1-0.5*(erf((x-muB)/sqrt(2*sigmaB^2))+erf((x+muB)/sqrt(2*sigmaB^2))))+2*log(2*sqrt(pi))+p(3),'--b')

plot(xB,polyval(p,xB),'r','LineWidth',4);

plot(xB,log(1-normcdf(xB,muB,sigmaB)))
plot(sort(data),log((length(data):-1:1)/length(data)),'--k','LineWidth',4);

% plot(x,log(0.5*(1-erf((x-muB)/sqrt(2*sigmaB^2))))+log(2*sqrt(pi)),'--b')

legend(['Folded-Normal        ';...
        'Folded-Normal O(x^2) ';...
        'Data                 '],'Location','SouthWest');
xlabel('B');
ylabel('log10(CDF^{-1})');
% axis([x(1) x(end) -3 0]);
% 
% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[-Inf,-1e-8],...
%                'Upper',[Inf,1e-8],...
%                'Startpoint',[1 1]);
% f = fittype('a*x+b','options',s);
% % [c2,gof]=fit(x1',y1',f);
% a=c2.a;
% plot(xB,a*xB.^2)
%%
%-------------------------------
figure;
axes('FontSize',24);
hold on;
grid on;
% scatter(log(abs(sinh(beta*SMFEnsemble/2))),log(2*abs(CurrentEstimateEnsemble)),ones(1e4,1)*0.1,'b');
% % scatter(log(abs(CurrentEstimateEnsemble)),-barrier,ones(1e4,1)*0.1);
% scatter(log(abs(sinh(beta*SMFEnsemble/2))),log(abs(CurrentEnsemble)),ones(1e4,1)*0.1,'r');
plot(log10(abs(CurrentEnsemble)),log10(abs(CurrentEstimateEnsemble)),'ob','MarkerSize',2);
% scatter(log(abs(CurrentEstimateEnsemble)),-barrier,ones(1e4,1)*0.1);
% plot(log(abs(CurrentEnsemble)),log(abs(sinh(beta*SMFEnsemble/2))),'.r','MarkerSize',2);

plot(log10(abs(CurrentEnsemble)),log10(abs(2*wAverageEnsemble.*sinh(beta*SMFEnsemble/2))),'or','MarkerSize',4);
plot(log10(abs(CurrentEnsemble)),log10(wAverageEnsemble),'og','MarkerSize',2);
plot(log10(abs(CurrentEnsemble)),log10(exp(-barrier/2)),'oc','MarkerSize',2);
x=linspace(min(log10(abs(CurrentEnsemble))),max(log10(abs(CurrentEnsemble))),100);
plot(x,x,'-k','LineWidth',4)
axis tight
xlabel('log10|I|');
ylabel('log10|y|');
h=legend(['y = w exp(-B/2)2sinh(SMF/2)';...
          'y = w 2sinh(SMF)           ';...
          'y = w                      ';...
          'y = exp(-B/2)              ';...
          'y = I                      '],'Location','Best');
        set(h, 'Color', 'none')

% print(gcf, '-depsc2', 'Scatter4');
%----------------------------
   

%%
data = (CurrentEnsemble);
% data = abs(data);
mu = mean(CurrentEnsemble);
sigma = std(CurrentEnsemble);
x=linspace(min(data),max(data),100);
p = polyfit(sort(data(data<0))',log((1:length(data(data<0)))/length(data(data<0))),1)

f1 = polyval(p,x(x<0));
f2 = polyval(p,x(x>0));
F = [exp(f1)/2,1-0.5*exp(-f2)];
%---------------------------------
figure;
axes('FontSize',24);
hold on;
grid on
plot(sort(data),(1:length(data)),'k','LineWidth',4);
plot(x,(F*numRealizations),'.-b','LineWidth',2)
plot(x,(normcdf(x,mu,sigma)*numRealizations),'--r','LineWidth',2)

legend(['Data       ';...
        'Exponential';...
        'Normal     '],'Location','SouthEast');
xlabel('I');
ylabel('CDF');
axis tight
% print(gcf, '-depsc2', 'CDFvsI');


data = abs(CurrentEnsemble);  abs(sinh(SMFEnsemble*beta/2));
p3 = polyfit(sort(data)',log((length(data):-1:1)/length(data)),1)
f3 = polyval(p3,x);
x=linspace(min(data),max(data),100);

%-----------------------
figure;
axes('FontSize',24);
hold on;grid on

plot(sort(data),(1:length(data)),'k','LineWidth',4);
plot(x,((1-exp(p3(1)*x))*numRealizations),'.-b','LineWidth',2);
plot(x,(numRealizations*0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2)))),'--r','LineWidth',2);
 
legend(['Data       ';...
        'Exponential';...
        'Normal     '],'Location','SouthEast');
xlabel('|I|');
ylabel('CDF');
axis([x(1) x(end) 0 (numRealizations)])
axis tight
% print(gcf, '-depsc2', 'CDFvs|I|');
%-------------------------------
p = polyfit(sort(data)',log10((length(data):-1:1)/length(data)),2);

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort((data)),log10(length(data):-1:1),'k','LineWidth',4);
plot(x,log10((exp(p3(1)*x))*numRealizations),'.-b','LineWidth',2);

plot(x,log10(numRealizations*(1-0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',2);
% plot(x,polyval(p,x)+log10(numRealizations),'g','LineWidth',4);

legend(['Data       ';...
        'Exponential';...
        'Normal     ';...
        'Parabola   '],'Location','SouthWest');
xlabel('|I_{sinai}|');
ylabel('log10(CDF^{-1})');
a=(sort(data));
% axis([x(1) a(numRealizations-10) 1 log10(numRealizations)])
% print(gcf, '-depsc2', 'log10iCDFvs|I|');
% print(gcf, '-depsc2', 'invCDF1');
% axis tight
%----------------------
%%


%%
% data = (wAverageEnsemble.*sinh(SMFEnsemble*beta/2));
data = CurrentEstimateEnsemble;
data = CurrentEnsemble./sinh(SMFEnsemble*beta/2);
% data = abs(data);
mu = mean(data);
sigma = std(data);
x=linspace(min(data),max(data),100);
p = polyfit(sort(data(data<0))',log((1:length(data(data<0)))/length(data(data<0))),1)

f1 = polyval(p,x(x<0));
f2 = polyval(p,x(x>0));

F = [exp(f1)/2,1-0.5*exp(-f2)];
%---------------------------------
figure;
axes('FontSize',24);
hold on;
grid on
plot(sort(data),(1:length(data)),'k','LineWidth',4);
plot(x,(F*numRealizations),'.-b','LineWidth',2)
plot(x,(normcdf(x,mu,sigma)*numRealizations),'--r','LineWidth',2)

legend(['Data       ';...
        'Exponential';...
        'Normal     '],'Location','SouthEast');
xlabel('w sinh(SMF)');
ylabel('CDF');
axis tight
% print(gcf, '-depsc2', 'CDFvsI');

% data=data(abs(data)>0.007);
mu = mean(data);
sigma = std(data);
data = abs(data);  
x=linspace(min(data),max(data),100);
p3 = polyfit(sort(data)',log((length(data):-1:1)/length(data)),1)
f3 = polyval(p3,x);

%-----------------------
figure;
axes('FontSize',24);
hold on;grid on

plot(sort(data),(1:length(data)),'k','LineWidth',4);
plot(x,((1-exp(p3(1)*x))*length(data)),'.-b','LineWidth',2);
plot(x,(length(data)*0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2)))),'--r','LineWidth',2);
 
legend(['Data       ';...
        'Exponential';...
        'Normal     '],'Location','SouthEast');
xlabel('|wsinh(SMF)|');
ylabel('CDF');
axis([x(1) x(end) 0 (length(data))])
axis tight
% print(gcf, '-depsc2', 'CDFvs|I|');
%-------------------------------
p = polyfit(sort(data)',log10((length(data):-1:1)/length(data)),2);

figure;
axes('FontSize',24);
hold on;grid on

plot(sort((data)),log10(length(data):-1:1),'k','LineWidth',4);
plot(x,log10((exp(p3(1)*x))*length(data)),'.-b','LineWidth',2);
plot(x,polyval(p,x)+log10(length(data)),'g','LineWidth',4);

plot(x,log10(length(data)*(1-0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',2);
legend(['Data       ';...
        'Exponential';...
        'Parabola   ';...
        'Normal     '],'Location','SouthWest');
xlabel('|I_{Est}|');
ylabel('log10(CDF^{-1})');
a=(sort(data));
axis([x(1) a(length(data)-10) 1 log10(length(data))])
%%
