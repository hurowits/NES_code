x=linspace(-20,-1,100);
figure;
axes('FontSize',24);
hold on;grid on
% plot(x,log10(normpdf(x,mean(I),std(I))),'r');
% plot(x,(1/2*((x-mean(I))/std(I))),'b');
plot(x,(1/2*((x-mean(I))/std(I)))-log10(2),'b','LineWidth',4);

% plot(x, ((1./x.^2).*(exp(-((1./x-mean(1./I))/sqrt(2)/std(1./I)).^2))/sqrt(2*pi*std(1./I)^2)),'k')
% plot(x, 2*((x-mean(I))/std(I)),'--k')

plot(x(1:end),log10(0.5*(1-erf((1./x  - mean(1./I))/(sqrt(2)*std(1./I))))),'g','LineWidth',4)
plot(sort(I),log10((1:length(I))/length(I)),'r','LineWidth',4)
% plot(x,log10((1./x).^2.*exp(-(1./x-mean(1./I)).^2/2/std(1./I)^2))/sqrt(2*pi*std(1./I)^2))
legend(['x^{1/2}           ';...
        'Inverse Normal CDF';...
        'Data              '],'Location','SouthEast');
xlabel('log10(|I|)');
ylabel('log10(CDF)')
% print(gcf, '-depsc2', 'Sinai_CDF_I');

%%
x=linspace(min(barrier),max(barrier),100);

figure;
axes('FontSize',24);
hold on;grid on

plot(x(1:end),(0.5*(1-erf((1./x  - mean(1./barrier))/(sqrt(2)*std(1./barrier))))),'--k','LineWidth',4)
plot(x, (0.5*(1+erf((log(x)-mean(log(barrier)))/sqrt(2*std(log(barrier))^2)))),'g','LineWidth',4);
plot(x(1:end),(normcdf(x,mean(barrier),std(barrier))),'b','LineWidth',4)
plot(sort(barrier),((1:length(barrier))/length(barrier)),'r','LineWidth',4);
% plot(x,log10((1./x).^2.*exp(-(1./x-mean(1./I)).^2/2/std(1./I)^2))/sqrt(2*pi*std(1./I)^2))
legend(['Inv-Normal';...
        'Log-Normal';...
        'Normal    ';...
        'Data      '],'Location','SouthEast');

xlabel('barrier');
ylabel('CDF')
print(gcf, '-depsc2', 'Barrier1');

figure;
axes('FontSize',24);
hold on;grid on

plot(x(1:end),log10(0.5*(1-erf((1./x  - mean(1./barrier))/(sqrt(2)*std(1./barrier))))),'--k','LineWidth',4)
plot(x, log10(0.5*(1+erf((log(x)-mean(log(barrier)))/sqrt(2*std(log(barrier))^2)))),'g','LineWidth',4);
plot(x(1:end),log10(normcdf(x,mean(barrier),std(barrier))),'b','LineWidth',4)
plot(sort(barrier),log10((1:length(barrier))/length(barrier)),'r','LineWidth',4);
% plot(x,log10((1./x).^2.*exp(-(1./x-mean(1./I)).^2/2/std(1./I)^2))/sqrt(2*pi*std(1./I)^2))
legend(['Inv-Normal';...
        'Log-Normal';...
        'Normal    ';...
        'Data      '],'Location','SouthEast');

xlabel('barrier');
ylabel('log10(CDF)')
print(gcf, '-depsc2', 'Barrier2');

figure;
axes('FontSize',24);
hold on;grid on

plot(x(1:end),log10(1-0.5*(1-erf((1./x  - mean(1./barrier))/(sqrt(2)*std(1./barrier))))),'--k','LineWidth',4)
plot(x, log10(1-0.5*(1+erf((log(x)-mean(log(barrier)))/sqrt(2*std(log(barrier))^2)))),'g','LineWidth',4);
plot(x(1:end),log10(1-normcdf(x,mean(barrier),std(barrier))),'b','LineWidth',4)
plot(sort(barrier),log10(1-(1:length(barrier))/length(barrier)),'r','LineWidth',4);
% plot(x,log10((1./x).^2.*exp(-(1./x-mean(1./I)).^2/2/std(1./I)^2))/sqrt(2*pi*std(1./I)^2))
legend(['Inv-Normal';...
        'Log-Normal';...
        'Normal    ';...
        'Data      '],'Location','SouthEast');
xlabel('barrier');
ylabel('log10(CDF^{-1})')
print(gcf, '-depsc2', 'Barrier3');