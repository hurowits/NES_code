clear all
close all
N_vec = 1e2;
beta = 0.1;
Delta_0 = 1;
w_beta=1;
Estrength_vec = logspace(-5,3,100);
I2=zeros(1,length(Estrength_vec));
Current=zeros(1,length(Estrength_vec));

% sigma_eff1b = zeros(1,N);
%%

for qN=1:length(N_vec)
    N = N_vec(qN);
    
    energy_vec = rand(1,N);
    delta_vec = Delta_0*[diff(energy_vec), energy_vec(N)-energy_vec(1)];
    % delta_vec = Delta_0*rand(1,N);
    sigma =6;
    
    %generate transition rates
    %     w_drv = genRand(sigma,randperm(N),N+1,'sigma');
    w_drv =exp(rand(1,N)*sigma-sigma/2);
    % w_drv = logspace(-sigma/2,sigma/2,N);
    % w_drv = w_drv(randperm(N));
    w_drv = w_drv*mean(1./w_drv);
    % w_drv = w_drv/mean(w_drv);
    % w_drv = rand(1,N)*sigma;
    for iE=1:length(Estrength_vec)
        h = waitbar(iE/length(Estrength_vec));
        
        Estrength=Estrength_vec(iE);
        w_pl_init = (Estrength^2*w_drv + 2*w_beta./(1+exp(beta*delta_vec)));
        w_mi_init = (Estrength^2*w_drv + 2*w_beta./(1+exp(-beta*delta_vec)));
        w_pl = w_pl_init;
        w_mi = w_mi_init;
        den=0;
        for t = 1:N
            for s=1:N
                den = den + prod(w_pl(s+1:N))*prod(w_mi(s:-1:2));
            end
            w_pl = circshift(w_pl,[0,1]);
            w_mi = circshift(w_mi,[0,1]);
        end
        
        I2(iE) = (prod(w_pl) - prod(w_mi))/den;
        I_0 = -beta/N*dot(w_drv,delta_vec)*Estrength_vec(1)^2;
        
        I = fminsearch(@(I) f_error(I,w_pl,w_mi),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
        I_0 = I;
        P(1) = 1;
        for q = 1:N
            P(q+1) = (w_pl(q)*P(q) - I) / w_mi(q);
        end
        
        P = P(1:end-1)./sum(P(1:end-1));
        Current(iE) = I;
        
    end
end