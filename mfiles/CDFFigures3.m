function [mu,sigma] =  CDFFigures3(data,fileName,xlab,flagPrint)
numRealizations = length(data);
mu = mean(data);
sigma = std(data);
x=linspace(min(data),max(data),100);

%----------------------------------------------%
x=linspace(min(data),max(data),100);
p3 = polyfit(sort(data)',log10((length(data):-1:1)/length(data)),1)
f3 = polyval(p3,x);

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort((data)),log10(length(data):-1:1),'-k','LineWidth',4);
plot(x,f3 + log10(numRealizations),'-b','LineWidth',4);
% plot(x,log10(numRealizations*(1-0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',4);
plot(x,log10(numRealizations*(1-normcdf(x,mu,sigma))),'--r','LineWidth',4);

legend(['Data       ';...
    'Exponential';...
%     'Parabola   ';...
    'Normal     '],'Location','NorthWest');
xlabel(xlab);
ylabel('log10(CDF^{-1})');
a=(sort(data));
axis([x(1) a(numRealizations) 1 log10(numRealizations)])
% print(gcf, '-depsc2', 'log10iCDFvs|I|');

if (flagPrint==1)
print(gcf, '-depsc2', fileName);
end
%----------------------------------------------%
x=linspace(min(data),max(data),100);
p3 = polyfit(sort(data)',log10((1:length(data))/length(data)),1)
f3 = polyval(p3,x);

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort((data)),log10(1:length(data)),'-k','LineWidth',4);
plot(x,f3 + log10(numRealizations),'-b','LineWidth',4);
% plot(x,log10(numRealizations*(0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',4);
plot(x,log10(numRealizations*(normcdf(x,mu,sigma))),'--r','LineWidth',4);

legend(['Data       ';...
    'Exponential';...
%     'Parabola   ';...
    'Normal     '],'Location','NorthWest');
xlabel(xlab);
ylabel('log10(CDF)');
a=(sort(data));
axis([x(1) a(numRealizations) 1 log10(numRealizations)])
% print(gcf, '-depsc2', 'log10iCDFvs|I|');

if (flagPrint==1)
print(gcf, '-depsc2', fileName);
end