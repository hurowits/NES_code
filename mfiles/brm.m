clear all

cols = 'rgbckym';c=0;
M=200;

beta=0.2;
w_beta=1;
a_0=1;L=a_0*M;
Delta=10;

s0 = (linspace(0,1,M)*2-1);
% s0 = randn(1,M);
randInd1 = randperm(M);
randInd2 = randperm(M);
% tau_vec=linspace(-1,1,100);       %Scaled driving intensity
nu_vec = logspace(-10,1,1000);
% nu_vec=linspace(0,100,1000);
% nu_vec = linspace(0,3e-6,1000);
% nu_vec=linspace(0,.01,100);
sigma_vec=5;
Energy = randn(1,M)*Delta;
% Energy = linspace(-Delta/2,Delta/2,M);
% Energy=Energy(randInd1);

Delta_n = Energy(2:M)-Energy(1:M-1);
Delta_n = [Delta_n, Energy(1)-Energy(M)];


w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:M-1)-Energy(2:M)))); %w_plus
w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(M)-Energy(1))))];

w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:M-1)-Energy(2:M)))); %w_minus
w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(M)-Energy(1))))];
SMF=zeros(1,length(nu_vec));
velocity=zeros(1,length(nu_vec));
Diffusion=zeros(1,length(nu_vec));
% for iSig = 1:length(sigma_vec);
%     sigma = sigma_vec(iSig);
sigma=sigma_vec(1);
G = 10.^(sigma*s0(randInd1));
G = G*mean(1./G);
% G=[1,0,0];
% nu_vec = 1/max(G)*10.^(sigma*tau_vec);
%%
% nu_vec = linspace(0,.1,1000);
velocity=zeros(1,length(nu_vec));
Diffusion=zeros(1,length(nu_vec));
SMF=zeros(1,length(nu_vec));
v_drda=zeros(1,length(nu_vec));
D_drda=zeros(1,length(nu_vec));
mu1=zeros(1,length(nu_vec));
mu2=zeros(1,length(nu_vec));
for q=1:length(nu_vec)
    nu = nu_vec(q);
    
    
    w_drv = nu*G;
    w_m = (w_drv + w_bath_ccw);
    w_p = (w_drv + w_bath_cw);
    
    
    
    
    
    SMF(q) = sum(log(w_p./w_m));
    
    s_eff(q) = log(max(w_m./w_p)) - log(min(w_m./w_p));
    
    [velocity(q),Diffusion(q)] = VD(w_p,w_m);
    
%     if (SMF(q)>0)
%         mu1(q) = mean(w_m./w_p);
%         mu2(q) = mean((w_m./w_p).^2);
%         if (mu1(q)<1)
%             v_drda(q) = (1- mean(w_m./w_p))./(mean(1./w_p));
%         else
%             %             v_drda(q) = 0;
%         end
%         D_drda(q) = 0.5*(mean(1./w_p)).^(-3)*...
%             (1-(mean(w_m./w_p)).^2)./...
%             (1-mean((w_m./w_p).^2))*...
%             (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
%             2*mean(w_m./w_p.^2)*mean(1./w_p));
%         
%         %         mu1(q) = mean(w_m./w_p);
%         %         mu2(q) = mean((w_m./w_p).^2);
%         %
%     else
%         mu1(q) = mean(w_p./w_m);
%         mu2(q) = mean((w_p./w_m).^2);
%         if(mu1(q)<1)
%             v_drda(q) = -(1- mean(w_p./w_m))./(mean(1./w_m));
%         else
%             %             v_drda(q)=0;
%         end
%         D_drda(q) = 0.5*(mean(1./w_m)).^(-3)*...
%             (1-(mean(w_p./w_m)).^2)./...
%             (1-mean((w_p./w_m).^2))*...
%             (mean(1./w_m.^2)*(1-mean(w_p./w_m))+...
%             2*mean(w_p./w_m.^2)*mean(1./w_m));
%         %
%         %         mu1(q) = mean(w_p./w_m);
%         %         mu2(q) = mean((w_p./w_m).^2);
%         
%         %         v_drda(q,iSig) = (1- mean(w_m./w_p))./(mean(1./w_p));
%         %         D_drda(q,iSig) = 0.5*(mean(1./w_p)).^(-3)*...
%         %             (1-(mean(w_m./w_p))f.^2)./...
%         %             (1-mean((w_m./w_p).^2))*...
%         %             (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
%         %             2*mean(w_m./w_p.^2)*mean(1./w_p));
%         %
%         
%     end
%     if (D_drda(q)<0)
%         D_drda(q)=0;
%     end
end %for q
a_inf = mean((1./w_p).^2)/mean(1./w_p)^2;
% end %for sigma
%         f(iR,:)=-(real(velocity(:,2))./real(Diffusion(:,2)))/2;
%     f(iR,:)=-(real(velocity(:,1))./real(Diffusion(:,1)))/2;

% end
%%
%
% figure;
% axes('FontSize',24);
% grid on
% hold on;
% plot(nu_vec, (real(velocity)),nu_vec,-SMF/M^2,'LineWidth',4);
% xlabel('\nu [scaled]');
% ylabel('V');
% axtype(1)
% axis([1e-7 5 0 10])
% print(gcf, '-depsc2', '/users/physics/hurowits/MyWebSite/PROJ/NES/Figs/DvsS_500.eps');
%%
figure;
axes('FontSize',24);
grid on
hold on;
plot(nu_vec, (v_drda./D_drda),'b','LineWidth',4);
plot(nu_vec, (velocity./Diffusion)/M,'g','LineWidth',4);

plot(nu_vec, SMF/M,'--r','LineWidth',4);
% plot(nu_vec,sqrt(sigma/2)*ones(size(SMF)),'--k','LineWidth',2)
% plot(nu_vec,-sqrt(sigma/2)*ones(size(SMF)),'--k','LineWidth',2)
xlabel('\nu');
ylabel('v/D');
legend(['derrida ';,...
    'numerics';...
    's       '],'Location','SouthEast');
% axis([nu_vec(1) nu_vec(end) -2e-3 1e-2])

axtype(1);
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/BRM/Figs/vD_derrida.eps')

%%

figure;
axes('FontSize',24);
grid on
hold on;
plot(nu_vec,v_drda,'b','LineWidth',4)

plot(nu_vec,-velocity*M,'--g','LineWidth',4);

% plot(nu_vec,SMF/M,'r','LineWidth',4);
axis([nu_vec(1) nu_vec(end) -0.02 0.06])

axtype(1)
xlabel('\nu')
ylabel('v');
legend(['derrida ';...
        'numerics'],'Location','NorthWest');
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/BRM/Figs/v_derrida.eps')
%%
figure;
axes('FontSize',24);
grid on
hold on;
plot(nu_vec,D_drda,'b','LineWidth',4)
plot(nu_vec,Diffusion*M^2,'--g','LineWidth',4);
% plot(nu_vec,abs(SMF)/M,'r','LineWidth',4);
xlabel('\nu');
ylabel('D')
axis([nu_vec(1) nu_vec(end) 0 50])
axtype(1);
legend(['derrida ';...
    'numerics'],'Location','NorthWest');

%     print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/BRM/Figs/D_derrida.eps')

%%
figure;
axes('FontSize',24);
grid on
hold on;
plot(nu_vec,mu1,'--b','LineWidth',4);
plot(nu_vec,mu2,'g','LineWidth',4);
% plot(nu_vec,SMF/M,'r','LineWidth',4)
plot(nu_vec,ones(size(nu_vec)),'--k');
xlabel('\nu')
legend(['\mu = 1';...
    '\mu = 2']);
axtype(1);
axis([nu_vec(1) nu_vec(end) 0.98 1.04]);
%     print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/BRM/Figs/derrida_mu.eps')
%%
%
% figure;
% grid on
% hold on;
% plot(nu_vec,v_drda,'b','LineWidth',4)
% plot(nu_vec,-velocity*M,'--b','LineWidth',4)
% plot(nu_vec,(D_drda/M^2),'r','LineWidth',4)
% plot(nu_vec,(Diffusion),'--r','LineWidth',4)
% plot(nu_vec,mu1,'-.b','LineWidth',4);
% plot(nu_vec,mu2,'-.r','LineWidth',4);
% axtype(1)