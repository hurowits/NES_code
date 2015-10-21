Nmax=1000;
numRealizations=1;
sigma = 6;
Estrength_vec = logspace(-sigma/2,sigma/2,1000);
Delta=1;
SMF_d = zeros(1,length(Estrength_vec));
SMF_g = zeros(1,length(Estrength_vec));
SMF = zeros(1,length(Estrength_vec));
zerosInSMF = zeros(1,Nmax);
Energy = randn(1,Nmax)*Delta;
for N = Nmax
    for iR = 1:numRealizations
        %     wb=waitbar(iR/numRealizations);
        wb=waitbar(N/Nmax);
        Energy = randn(1,N)*Delta;
        Delta_n = zeros(1,N);
        Delta_n(1:N-1) = Energy(2:N)-Energy(1:N-1);
        Delta_n(N) =  Energy(1)-Energy(N);
        % G =10.^(rand(1,N)*sigma - sigma/2);
        G =10.^(linspace(-sigma/2,sigma/2,N));
        G=G(randperm(N));
        harmonicMeanG= mean(1./G)^-1;
        % G = G*mean(1./G);
        
        for iE = 1:length(Estrength_vec)
            SMF(iE) = sum(Delta_n./(1+Estrength_vec(iE)^2*G));
            SMF_d(iE) = sum (Delta_n(G<1/Estrength_vec(iE)^2));
            SMF_g(iE) = 0.5*sum(Delta_n.*erfc(log10(G) +log10(Estrength_vec(iE)^2)));
            
        end
        zerosInSMF(N) = sum(abs(diff(sign(SMF))))/2;
    end
end
%%
figure;
axes('FontSize',24);
plot(1:Nmax,zerosInSMF,'.');
hold on;
plot(1:Nmax, sqrt(pi*sigma)*ones(1,Nmax),'r-')
xlabel('Length of chain');
ylabel('Number of sign changes');

figure;
axes('FontSize',24);
plot(log10(Estrength_vec.^2),((SMF_d)),'b','LineWidth',3)
hold on
plot(log10(Estrength_vec.^2),((SMF)),'r','LineWidth',3)
grid on
plot(log10(Estrength_vec.^2),((SMF_g)),'g','LineWidth',3)

xlabel('log_{10}\epsilon^2')
ylabel('SMF');
    axis([-sigma/2 sigma/2 min([SMF';SMF_d']) max([SMF';SMF_d'])])

legend(['Discrete  ';...
    %         'Continuous';...
    'erfc      ']);
% print(gcf, '-depsc2', 'SMF_RW.eps');

%%
prefactor = zeros(length(Estrength_vec),length(G));
prefactor_d = zeros(length(Estrength_vec),length(G));
legendText=[];
Gsorted = sort(G);
iGG=0;
for iG = 1:300:1000
    iGG=iGG+1;
    prefactor(:,iGG) = 1./(1+Gsorted(iG)*Estrength_vec.^2);
    prefactor_d(:,iGG) = Gsorted(iG)<1./Estrength_vec.^2;
    %     prefactor_g(:,iGG) = erfc(Gsorted(iG)*Estrength_vec.^2/sqrt(2));
    %     prefactor_g(:,iGG) =2./(1+exp(Gsorted(iG)*Estrength_vec.^2));
    prefactor_g(:,iGG) =0.5*erfc((log10(Estrength_vec.^2)+log10(Gsorted(iG))) );
    
    
    %     =[legendText;num2str(G(iG))];
end
legendText = num2str(Gsorted(1:300:1000)','%10.2e');

figure;
axes('FontSize',24);
% hold on
semilogx((Estrength_vec.^2), prefactor,'LineWidth',4);
hold on
plot((Estrength_vec.^2),prefactor_d,'-.','LineWidth',4);
plot((Estrength_vec.^2),prefactor_g,'--','LineWidth',2);
plot(1./(Gsorted(1:300:1000)),0.5*ones(1,length(Gsorted(1:300:1000))),'dg','LineWidth',5)
grid on
axis([10^(-sigma) 10^(sigma) 0 1])
xlabel('\epsilon^2');
ylabel('[1+g\epsilon^2]^{-1}')
legend(legendText);
% print(gcf, '-depsc2', 'prefactor');
%%
N=5;
r=1:N;
p = 2^(-2*N)* factorial(2*N+1)./factorial(N-r)./factorial(N+r+1);
figure;plot(r,(p),'.')