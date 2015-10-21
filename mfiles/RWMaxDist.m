n=30;

P1=zeros(n,1);
P2=zeros(n,1);
for k=1:n
    for i=1:(n/k+1)/2
        P1(k) = P1(k) + factorial(2*n)/(factorial(2*n-(n-(2*i-1)*k))*factorial(n-(2*i-1)*k));
    
    end
end

for k=1:n
    for i=1:(n/k)/2
     P2(k) = P2(k) + factorial(2*n)/(factorial(2*n-(n-2*i*k))*factorial(n-2*i*k));
    end
end
    
P = 1- 2*(P1-P2)/(factorial(2*n)/(factorial(n)^2));
figure;
plot(1:n,log10(P),'-')