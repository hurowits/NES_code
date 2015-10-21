clear all
M=100;

sigma_vec = linspace(0,10,200);
SMF_vec = linspace(0,max(sigma_vec)*M,100);
Delta = 5
beta=10
w_beta=1
% nu_vec = 
Energy = randn(1,M)*Delta;
Energy = linspace(-Delta,Delta,M)
Delta_n = Energy(2:M)-Energy(1:M-1);
Delta_n = [Delta_n, Energy(1)-Energy(M)];


w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:M-1)-Energy(2:M)))); %w_plus
w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(M)-Energy(1))))];

w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:M-1)-Energy(2:M)))); %w_minus
w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(M)-Energy(1))))];

G=1;
s0 = (linspace(0,1,M)*2-1);
randInd1 = randperm(M);

for iSigma =1:length(sigma_vec)
    sigma = sigma_vec(iSigma);
    
    for iSMF = 1:length(SMF_vec)
        SMF = SMF_vec(iSMF);
        
        
        s_p = sigma*s0(randInd1)+SMF/M*ones(1,M);
        
        w_p = G.*exp(s_p/2);
        w_m = G.*exp(-s_p/2);
        mu1(iSigma,iSMF) = mean(w_m./w_p);
        mu2(iSigma,iSMF) = mean((w_m./w_p).^2);
    end
end

G = 10.^(sigma*s0(randInd1));
G = G*mean(1./G);

    nu_vec = logspace(-10,10,1000);

for iSigma =1:length(sigma_vec);
    sigma = sigma_vec(iSigma);

    for q = 1:length(nu_vec)
      nu = nu_vec(q);
    SMF(iSigma,q) = sum(log(w_p./w_m));
    
    w_drv = nu*G;
    w_m = (w_drv + w_bath_ccw);
    w_p = (w_drv + w_bath_cw);
    s_eff(iSigma,q)=log(max(w_m./w_p)./min(w_m./w_p));
      if(SMF(iSigma,q)>0)
        mu1_brm(iSigma,q) = mean(w_m./w_p);
        mu2_brm(iSigma,q) = mean((w_m./w_p).^2);
      else
           mu1_brm(iSigma,q) = mean(w_p./w_m);
        mu2_brm(iSigma,q) = mean((w_p./w_m).^2);
      end
    end
end
%%
figure;imagesc(SMF_vec/M,sigma_vec,(mu1)>1)
xlabel('s');ylabel('\sigma');
title(['\mu_1, M=',num2str(M)]);


figure;imagesc(SMF_vec/M,sigma_vec,(mu2)>1)
xlabel('s');ylabel('\sigma');
title(['\mu_2, M=',num2str(M)]);

%%
figure;imagesc(log10(nu_vec),sigma_vec,(SMF))

figure;imagesc(log10(nu_vec),sigma_vec,mu1_brm>1)

figure;imagesc(log10(nu_vec),sigma_vec,mu2_brm>1)
