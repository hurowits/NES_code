%% initialization

clear all
close all
methodRG =1;
%1 = Saar, 2 = Monthous
N_vec = 1e6;[1e3,1e4,1e5,1e6,1e7]; %number of sites
nBins = 10000;
beta = 1;
Delta_0 = 1;
w_beta=1;
alpha=0;

sigma_eff1 = zeros(length(N_vec),max(N_vec));
% sigma_eff1b = zeros(1,N);
%%

for qN=1:length(N_vec)
    N = N_vec(qN);
    h = waitbar(0,['renormalizing ',num2str(N), 'site ring']);

    energy_vec = rand(1,N);
    delta_vec = Delta_0*[diff(energy_vec), energy_vec(N)-energy_vec(1)];
    % delta_vec = Delta_0*rand(1,N);
    sigma =3;
    
    %generate transition rates
%     w_drv = genRand(sigma,randperm(N),N+1,'sigma');
    w_drv =exp(rand(1,N)*sigma-sigma/2);
    % w_drv = logspace(-sigma/2,sigma/2,N);
    % w_drv = w_drv(randperm(N));
    w_drv = w_drv*mean(1./w_drv);
    % w_drv = w_drv/mean(w_drv);
    % w_drv = rand(1,N)*sigma;
    w_pl_init = log(w_drv + 2*w_beta./(1+exp(beta*delta_vec)));
    w_mi_init = log(w_drv + 2*w_beta./(1+exp(-beta*delta_vec)));
    coupling = log(exp(w_pl_init)+exp(w_mi_init));
    escapeRate = log(exp(w_pl_init) + exp(w_mi_init([N,1:N-1])));
    gamma1 = zeros(1,N);
    gamma2 = zeros(1,N);
    span1 = zeros(1,N);
    span1b = zeros(1,N);
    w_pl = w_pl_init;
    w_mi = w_mi_init;
    
   
    for r = 1:(N)
        
        waitbar(r/N);
%         tic;
        [maxVal,I]=max(coupling);
        I = I-1;
        sigma_eff1(qN,r) = std((coupling));
        %             sigma_eff1b(qN,r) = std(exp(coupling));
        gamma1(r) = maxVal;
        
        Imi = mod(I-1,N-r+1)+1;
        Ipl = mod(I+1,N-r+1)+1;
       
        w_pl(Ipl) = w_pl(I+1)+w_pl(Ipl)-maxVal;
        w_mi(Imi) = w_mi(I+1)+w_mi(Ipl)-maxVal;
        
        w_pl(I+1)   = [];
        w_mi(I+1) = [];
        coupling = log(exp(w_pl)+exp(w_mi));
%         toc;
        
    end
    close(h);
%     save ('RGData');
%     figure;plot(1:N,sigma_eff1(qN,1:N))

end


figure(1);
axes('FontSize',20);
hold on; 
plot((1:N_vec(1))/N_vec(1),sigma_eff1(1,1:N_vec(1)),'r','LineWidth',2);
plot((1:N_vec(2))/N_vec(2),sigma_eff1(2,1:N_vec(2)),'g','LineWidth',2);
plot((1:N_vec(3))/N_vec(3),sigma_eff1(3,1:N_vec(3)),'b','LineWidth',2);
plot((1:N_vec(4))/N_vec(4),sigma_eff1(4,1:N_vec(4)),'c','LineWidth',2);
plot((1:N_vec(5))/N_vec(5),sigma_eff1(5,1:N_vec(5)),'k','LineWidth',2);

legend(['N = 10^3';...
        'N = 10^4';...
        'N = 10^5';...
        'N = 10^6';...
        'N = 10^7'])
xlabel('\alpha = M/N');
ylabel('\sigma');
hold off
print(gcf, '-depsc2', ['RG_width_N.eps']);


figure(4);
axes('FontSize',20);
% plot(1:r,sigma_eff);
plot(1:r,(sigma_eff1),'r','LineWidth',2);hold on
plot(1:r,(sigma_eff2),'b','LineWidth',2);hold off
legend(['method 1';'method 2'],'Location','Best')
xlabel('# RG steps');
% xlabel('RG steps');
ylabel('std [ log(w) ]');
% print(gcf, '-depsc2', ['RG_spreading1.eps']);
