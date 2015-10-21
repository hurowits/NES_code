%calculate pdf of maximum value of random walk 
N=1000
k=linspace(0,N-1,1000);

p = log((N+k)./(N-k))./(1-k.^2/N^2).^N; 

lnp = log(log((N+k)./(N-k))) - gammaln(N+k)-gammaln(N-k);

p=exp(lnp/100);
C1=sum(p);
p=(p/C1).^100
p=p./sum(p);
figure;
plot(k,p)
