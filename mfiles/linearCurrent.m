N=1000;
num = 10000;
sigma = 5;

g_n = 10.^(rand(num,N)*sigma-sigma/2);

w= 1./(1/N * sum(1./g_n,2));
w_perm = w(randperm(num)');
harmonicMeanG = mean(1./g_n,2).^-1;
g_n = g_n./repmat(harmonicMeanG,1,N);

Energy = randn(1,N);

Delta_n = Energy(2:N)-Energy(1:N-1);
Delta_n = [Delta_n, Energy(1)-Energy(N)];
Delta_n = repmat(Delta_n,num,1);
Delta0 = sum(Delta_n.*g_n,2);
DeltaInf =  sum(Delta_n./g_n,2);
% I= Delta0./sum(1./g_n,2)/N;
I= Delta0.*w*beta/N;
I_perm = Delta0.*w_perm/N;
I_inf = DeltaInf*beta/N;
Imin = min(I);
Imax = max(I);
logImin = min(log10(abs(I)));
logImax = max(log10(abs(I)));

I_inf_min = min(log10(abs(I_inf)));
I_inf_max = max(log10(abs(I_inf)));


x_I = linspace(Imin,Imax,100);
x_logI = linspace(logImin,logImax,100);

x_I_inf = linspace(I_inf_min,I_inf_max,100);

x_Delta0 = linspace(min(Delta0),max(Delta0));
%%
figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(DeltaInf),((1:length(DeltaInf))/length(DeltaInf)),'r','LineWidth',4);
plot(linspace(min((DeltaInf)),max((DeltaInf)),100),(normcdf(linspace(min(DeltaInf),max(DeltaInf),100),mean(DeltaInf),std(DeltaInf))),'LineWidth',2);
legend(['Data        ';...
        'normal dist.';...
        'x^{1/2}     ']);

xlabel('Delta^{inf}');
ylabel('CDF');
%  print(gcf, '-depsc2', 'DeltaInf1');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(DeltaInf),log10((1:length(DeltaInf))/length(DeltaInf)),'r','LineWidth',4);
plot(linspace(min((DeltaInf)),max((DeltaInf)),100),log10(normcdf(linspace(min(DeltaInf),max(DeltaInf),100),mean(DeltaInf),std(DeltaInf))),'LineWidth',2);
legend(['Data        ';...
        'normal dist.';...
        'x^{1/2}     ']);

xlabel('Delta^{inf}');
ylabel('CDF');
% print(gcf, '-depsc2', 'DeltaInf2');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(DeltaInf),log10(1-(1:length(DeltaInf))/length(DeltaInf)),'r','LineWidth',4);
plot(linspace(min((DeltaInf)),max((DeltaInf)),100),log10(1-normcdf(linspace(min(DeltaInf),max(DeltaInf),100),mean(DeltaInf),std(DeltaInf))),'LineWidth',2);
legend(['Data        ';...
        'normal dist.';...
        'x^{1/2}     ']);

xlabel('Delta^{inf}');
ylabel('CDF');
% print(gcf, '-depsc2', 'DeltaInf3');
%% %%
figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(I),((1:length(I))/length(I)),'r','LineWidth',4);
% plot(sort(log10(abs(I_perm))),log10((1:length(I))/length(I)),'--k','LineWidth',4);
plot(linspace(min((I)),max((I)),100),(normcdf(linspace(min((I)),max((I)),100),mean(((I))),std(((I))))),'LineWidth',2);

legend(['Data        ';...
        'normal dist.';...
        'x^{1/2}     ']);

xlabel('I_0');
ylabel('CDF');
% print(gcf, '-depsc2', 'I0_1');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(I),log10((1:length(I))/length(I)),'r','LineWidth',4);
plot(x_I,log10(normcdf(x_I,mean(I),std(I))),'LineWidth',2);
% plot(x_I,  0.5*(x_I-mean(I)/std(I)),'--','LineWidth',2);

legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);

xlabel('I_0');
ylabel('log10(CDF)');
% print(gcf, '-depsc2', 'I0_2');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(I),log10(1-(1:length(I))/length(I)),'r','LineWidth',4);
plot(x_I,log10(1-normcdf(x_I,mean(I),std(I))),'LineWidth',2);
% plot(x_I, 1-0.5*((x_I-mean(log10(abs(I))))/std(log10(abs(I))))-log10(2),'--','LineWidth',2);

legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);

xlabel('I_0');
ylabel('log10(CDF^{-1})');
% print(gcf, '-depsc2', 'I0_3');



%%
figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(log10(abs(I))),((1:length(I))/length(I)),'r','LineWidth',4);
% plot(sort(log10(abs(I_perm))),log10((1:length(I))/length(I)),'--k','LineWidth',4);
plot(x_logI,normcdf(x_logI,mean(log10(abs(I))),std(log10(abs(I)))),'LineWidth',2);


legend(['Data        ';...
        'normal dist.';...
        'x^{1/2}     ']);

xlabel('I_0');
ylabel('CDF');
% print(gcf, '-depsc2', 'log|I0|1');

%%
figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(log10(abs(I))),log10((1:length(I))/length(I)),'r','LineWidth',4);
plot(x_logI,log10(normcdf(x_logI,mean(log10(abs(I))),std(log10(abs(I))))),'LineWidth',2);
plot(x_logI,  0.5*((x_logI-mean(log10(abs(I))))/std(log10(abs(I))))-log10(2),'--','LineWidth',2);
% plot(x_logI,log10(1-0.5*(1-erf((1./x_logI  - mean(1./log10(abs(I))))/(sqrt(2)*std(1./log10(abs(I))))))),'g:','LineWidth',2)
plot(x_logI,((0.5*(1-erf((1./x_logI  - mean(1./log10(abs(I))))/(sqrt(2)*std(1./log10(abs(I)))))))),'g:','LineWidth',2)

legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);

xlabel('log10|I_0|');
ylabel('log10(CDF)');
% print(gcf, '-depsc2', 'log|I0|2');



figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(log10(abs(I))),log10(1-(1:length(I))/length(I)),'r','LineWidth',4);
plot(x_logI,log10(1-normcdf(x_logI,mean(log10(abs(I))),std(log10(abs(I))))),'LineWidth',2);
% plot(x_I, 1-0.5*((x_I-mean(log10(abs(I))))/std(log10(abs(I))))-log10(2),'--','LineWidth',2);

legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);

xlabel('log10|I_0|');
ylabel('log10(CDF^{-1})');
% print(gcf, '-depsc2', 'log|I0|3');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(log10(abs(I_inf))),log10((1:length(I_inf))/length(I_inf)),'r','LineWidth',4);
plot(x_I_inf,log10(normcdf(x_I_inf,mean(log10(abs(I_inf))),std(log10(abs(I_inf))))),'LineWidth',2);
plot(x_I_inf, 0.5*((x_I_inf-mean(log10(abs(I_inf))))/std(log10(abs(I_inf))))-log10(2),'--','LineWidth',2);
legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);
    
xlabel('log10|I_{sat}|');
ylabel('log10(CDF)');



figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(log10(abs(I_inf))),log10(1-(1:length(I_inf))/length(I_inf)),'r','LineWidth',4);
plot(x_I_inf,log10(1-normcdf(x_I_inf,mean(log10(abs(I_inf))),std(log10(abs(I_inf))))),'LineWidth',2);
legend(['Data   ';...
        'Normal ';...
        'x^{1/2}']);
    
xlabel('log10|I_{sat}|');
ylabel('log10(CDF^{-1})');
%% Delta0
figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(Delta0),((1:length(I))/length(I)),'r','LineWidth',2)
plot(x_Delta0,(normcdf(x_Delta0,mean(Delta0),std(Delta0))),'LineWidth',2);
legend(['Data  ';...
        'Normal']);
xlabel('Delta^{(0)}');
ylabel('(CDF)');
% print(gcf, '-depsc2', 'Delta0_1');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(Delta0),log10((1:length(I))/length(I)),'r','LineWidth',2)
plot(x_Delta0,log10(normcdf(x_Delta0,mean(Delta0),std(Delta0))),'LineWidth',2);
legend(['Data  ';...
        'Normal']);
xlabel('Delta^{(0)}');
ylabel('log10(CDF)');
% print(gcf, '-depsc2', 'Delta0_2');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(Delta0),log10(1-(1:length(I))/length(I)),'r','LineWidth',2)
plot(x_Delta0,log10(1-normcdf(x_Delta0,mean(Delta0),std(Delta0))),'LineWidth',2);
legend(['Data      ';...
    'log-normal']);
xlabel('Delta^{(0)}');
ylabel('log10(CDF^{-1})');
% print(gcf, '-depsc2', 'Delta0_3');

%% w


wmin = min(w);
wmax = max(w);

x_w = linspace(wmin,wmax,100);
mu=mean(1./w);
sig = std(1./w);

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(w),((1:length(w))/length(w)),'r','LineWidth',2)
plot(x_w,(normcdf(x_w,mean(w),std(w))),'LineWidth',2);
plot(x_w,(0.5*(1-erf((1./x_w - mean(1./w))/(sqrt(2)*std(1./w))))),'--','LineWidth',2)
% plot(x_G, (0.5*(1+erf(mu/sig/sqrt(2)) + 1./x_G*sqrt(2/pi))),'g')

legend(['Data          ';...
        'log-normal    ';...
        'Inverse normal'],'Location','SouthEast');
xlabel('w^{\epsilon}');
ylabel('(CDF)');
% print(gcf, '-depsc2', 'w1');


figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(w),log10((1:length(w))/length(w)),'r','LineWidth',2)
plot(x_w,log10(normcdf(x_w,mean(w),std(w))),'LineWidth',2);
plot(x_w,log10(0.5*(1-erf((1./x_w - mean(1./w))/(sqrt(2)*std(1./w))))),'--','LineWidth',2)
% plot(x_G, log10(0.5*(1+erf(mu/sig/sqrt(2)) + 1./x_G * exp(-mu^2/2/sig^2)*sqrt(2/pi/sig^2))))

legend(['Data          ';...
        'log-normal    ';...
        'Inverse normal'],'Location','SouthEast');
xlabel('w^{\epsilon}');
ylabel('log10(CDF)');
% print(gcf, '-depsc2', 'w2');

figure;
axes('FontSize',24);
hold on;
grid on;
plot(sort(w),log10(1-(1:length(I))/length(I)),'r','LineWidth',2)
plot(x_w,log10(1-normcdf(x_w,mean(w),std(w))),'LineWidth',2);
plot(x_w,log10(1-0.5*(1-erf((1./x_w - mean(1./w))/(sqrt(2)*std(1./w))))),'--','LineWidth',2)
% plot(x_G, log10(1-0.5*(1+erf(mu/sig/sqrt(2)) + 1./x_G * exp(-mu^2/2/sig^2)*sqrt(2/pi/sig^2))),'g')


legend(['Data          ';...
        'log-normal    ';...
        'Inverse normal'],'Location','SouthWest');
xlabel('w^{\epsilon}');
ylabel('log10(CDF^{-1})');
% print(gcf, '-depsc2', 'w3');
