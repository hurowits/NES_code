xlabels=['I_{Sinai}     ';...
         'I_{linear}    ';...
         'I_{saturation}']
l=0;
for data = [CurrentEnsemble I I_inf]
    l=l+1;
% data = CurrentEnsemble;
% data = I_inf
% data = I;
mu = mean(data);
sigma = std(data);
data = abs(data);

p3 = polyfit(sort(data)',log((length(data):-1:1)/length(data))+log10(numRealizations),1)
f3 = polyval(p3,x);
x=linspace(min(data),max(data),100);

%-------------------------------

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort((data)),log10(length(data):-1:1),'k','LineWidth',4);
plot(x,log10((exp(p3(1)*x))*numRealizations),'.-b','LineWidth',4);
plot(x,log10(numRealizations*(1-0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',4);
legend(['Data         ';...
        'Exponential  ';...
        'Folded Normal'],'Location','SouthWest');
xlabel(['|',xlabels(l,:),'|']);
ylabel('log10(CDF^{-1})');
a=(sort(data));
axis([x(1) a(numRealizations-10) 1 log10(numRealizations)]);
% print(gcf, '-depsc2', ['invCDF',num2str(l)]);

end