tt = linspace(0,t_end,100);
figure;
hold on
plot(t_vec,site_t);
plot(tt,stdX*sqrt(tt)/N,'k','LineWidth',4)
plot(tt,-stdX*sqrt(tt)/N,'k','LineWidth',4);

tt=(1:ii-1)*dt;
figure;
hold on
plot(tt,stdX.^2/2);
plot(tt,tt,'g','LineWidth',4)

%
figure;
x = min(site_t(mask)):max(site_t(mask));
hist(site_t(mask),x);
hold on
plot(x,normpdf(x,mean(site_t(mask)),std(site_t(mask)))*sum(mask(:)),'r')
% plot(x,normpdf(x,mean(site_t(mask)),1/(1+1/100)/2)*sum(mask(:)),'r')