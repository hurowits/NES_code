% nBins = 100;
% mD=mean(log(D));
% mB=mean(log(B));
% mA=mean(log(A));
% mAB=mean(log(abs(B-A)));
% 
% [hA,xA] = hist(log(A),nBins);
% [hB,xB] = hist(log(B),nBins);
% [hAB,xAB] = hist(log(abs(A-B)),nBins);
% [hD,xD] = hist(log(D),nBins);
% 
% 
% 
% figure;plot(xA-(mA),hA,xB-(mB),hB,xD-(mD),hD,xAB-(mAB),hAB);
% legend(['A  ';'B  ';'D  ';'A-B'])
% % sum((A-mA).*(D-mD))./(sqrt(sum((A-mA).^2))*sqrt(sum((D-mD).^2)));
% % sum((B-mB).*(D-mD))./(sqrt(sum((B-mB).^2))*sqrt(sum((D-mD).^2)));
% % 
% % sum((A-B-mB-mA).*(D-mD))./(sqrt(sum((A-B-mB-mA).^2))*sqrt(sum((D-mD).^2)))
% 
% Current = (A-B)./D;

%%
nBins=50;
iE1 = 1;
a1 = std(Energy)*Estrength_vec(1)^2*beta/N;
f1 =  hist(Current(:,iE1)/a1,nBins);
f1=f1/sum(f1);
[hI,x1] = hist(Current(:,iE1),nBins);

iE2 = 30;
a2 = std(Energy)*Estrength_vec(iE2)^2*beta/N;
f2 =  hist(Current(:,iE2),nBins);
f2=f2/sum(f2);
[hI,x2] = hist(Current(:,iE2),nBins);

iE3 = 45;
a3 = std(Energy)*Estrength_vec(iE3)^2*beta/N;
f3 =  hist(Current(:,iE3),nBins);
f3=f3/sum(f3);
[hI,x3] = hist(Current(:,iE3),nBins);



x = linspace(min(x1),max(x1),nBins);
f = (x(2)-x(1))*exp(-(x).^2/ (2*std(Current(:,iE1))^2))/sqrt(2*pi*std(Current(:,iE1))^2);
% f = f/sum(f);
% 
% I0 =   std(Energy)*beta;
% % x = linspace(-10,10,nBins);
% f = (x(2)-x(1))*exp(-(x).^2/ (2*I0^2))/sqrt(2*pi*I0^2);
% f=f/sum(f);

figure;
axes('FontSize',20);
plot(x1/a1,f1,'d','MarkerSize',10);
grid
hold on;
plot(x2/a2,f2,'og','MarkerSize',10);
plot(x3/a3,f3,'.r','MarkerSize',15);
plot(x/a1,f,'k','LineWidth',2);
xlabel('I/I_0');
ylabel('P(I)');
l = legend(['\epsilon^2 < g_{max}^{-1}';...
        '\epsilon^2 < g_{max}^{-1}';...
        '\epsilon^2 > g_{max}^{-1}';...
        'Gaussian fit             ']);
    
%  print(gcf, '-depsc2', 'fvsI50.eps')


% figure;plot(x1/a1,f1,x2/a2,f2,x3/a3,f3, x1/a1,f)
%%
iE = 10;
[hI,xI] = hist(Current(:,iE),100);


hI=hI/sum(hI);

I_std_linear =  std(Energy)*Estrength_vec(iE)^2*beta
I_mean_linear = mean(G)*mean(Delta_n)*beta*Estrength_vec(iE)^2;
P_est_linear = exp(-(xI-I_mean_linear).^2/(2*I_std_linear^2));
P_est_linear = P_est_linear/sum(P_est_linear);

% 
% I_std_sat = sqrt(var(1./G)*2*N) * std(Delta_n)*Estrength_vec(iE)^2*beta/N
% I_mean_sat = mean(1./G)*mean(Delta_n)*beta*Estrength_vec(iE)^2
% P_est_sat = exp(-(xI-I_mean_sat).^2/(2*I_std_sat^2));
% P_est_sat = P_est_sat/sum(P_est_sat);

% mI = I_mean_linear;
mI = mean(Current(:,iE));
% sigI = I_std_linear;%
sigI=  std(Current(:,iE));
P_est = exp(-(xI-mI).^2/(2*sigI^2));
P_est = P_est/sum(P_est);

x_scaled = xI*sigI;
P_scaled =  sigI*exp(-(x_scaled-mI).^2/(2*sigI^2));
P_scaled = P_scaled/sum(P_scaled);

figure(1);
plot(xI,((hI)/sum(hI)),'.r','MarkerSize',10);hold on
plot(xI,((P_est)),'LineWidth',2);hold on
% plot(xI,((P_scaled)),'g','LineWidth',2);hold on

title('histogram of currents');
hold off

figure(2);
axes('FontSize',24);
hold on;
grid on;
plot(-(xI-mI).^2/(2*sigI^2),log(hI) ,'.','MarkerSize',10);hold on
plot(-(xI-mI).^2/(2*sigI^2) ,-(xI-mI).^2/(2*sigI^2) -  abs(max(log(hI))),'r');
xlabel('-(x-\mu)/2\sigma^2');
ylabel('ln[P(I)]');
hold off


figure(3);
axes('FontSize',24);
plot(xI,log(hI),'.','MarkerSize',10);hold on
plot(xI,(log(P_est)),'LineWidth',2)
xlabel('I');
ylabel('P(I)');
grid on;
hold off
