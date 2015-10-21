a = 1;b=1;
[s1,s2] = meshgrid( linspace(-20,20,100),linspace(-20,20,100));

% s=linspace(0,10,1000);
% s1 = zeros(1,1000);
% s2=s-s1;
sigma = 1- 4*a*b * cosh((s1+s2)/4).^2./(a*cosh(s1/2)+b*cosh(s2/2)).^2;
    rat2 = 4*tanh((s1+s2)/4)/(1+sigma * tanh((s1+s2)/4)^2);

figure;imagesc(linspace(-20,20,100),linspace(-20,20,100),rat2)
% 
% figure;surf(s1,s2,sigma)

% figure;plot(s,1./(coth(s/4).*(1+ sigma./coth(s/4).^2)/4))