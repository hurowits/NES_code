data = (CurrentEnsemble);
% data = abs(data);
mu = mean(CurrentEnsemble);
sigma = std(CurrentEnsemble);
x=linspace(min(data),max(data),100);

p = polyfit(sort(data(data<0))',log((1:length(data(data<0)))/length(data(data<0))),1)

f1 = polyval(p,x(x<0));
f2 = polyval(p,x(x>0));

F = [exp(f1)/2,1-0.5*exp(-f2)];

figure;
axes('FontSize',24);
hold on;grid on

plot(sort(data),(1:length(data)),'k','LineWidth',4);
plot(x,(F*numRealizations),'b','LineWidth',2)
plot(x,(normcdf(x,mu,sigma)*numRealizations),'g','LineWidth',2)

legend(['Data       ';...
        'f=exp(-|I|)';...
        'Normal     '],'Location','NorthWest');
xlabel('I');
ylabel('CDF');
% axis([x(1) x(end) 0 log10(numRealizations)])


data = abs(CurrentEnsemble);
p3 = polyfit(sort(data)',log((length(data):-1:1)/length(data)),1)
f3 = polyval(p3,x);
x=linspace(min(data),max(data),100);
figure;
axes('FontSize',24);
hold on;grid on

plot(sort((data)),(1:length(data)),'k','LineWidth',4);
plot(x,((1-exp(f3))*numRealizations),'b','LineWidth',2);
plot(x,(numRealizations*0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2)))),'g','LineWidth',2)
 
legend(['Data       ';...
        'exponential';...
        'Normal     '],'Location','NorthWest');
xlabel('|I|');
ylabel('CDF');
axis([x(1) x(end) 0 (numRealizations)])
