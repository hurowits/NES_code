N=1000
R= (randn(N,numRealizations));
R(R>0)=1;
R(R<0)=-1;
y = cumsum(R);
B = max(y)-min(y);
B1 = max(y)-mean(y);

x=linspace(0,1,numRealizations/100);

fB = 2*exp(-N*x.^2).*(sqrt(pi*N/2) .* exp(x.^2*N/2).*(x.^2*N-1).*erf(x*sqrt(N/2))+N*x);
fB=fB./sum(fB);
h = hist(B,x*N/sqrt(2));

figure;
plot(x,h/sum(h))
hold on
plot(x,fB,'r')
