% Script for calculating current and statistics for nearest neighbor jump
% process on a ring.

clear all

% -----------------------
% Running parameters
% -----------------------

numRealizations=1;
flagAnalytic=0;  %0 for numerical solution, 1 for analytic solution

% -----------------------
% Model parameters
% -----------------------

N=100; %number of sites

Delta = 1;                             %Energy interval
w_beta = 1;                            %Coupling to bath
beta = 0.1;                            %Inverse temperature
sigma_vec =3;                        %log-width of distribution


nu_vec= logspace(-8,5,100);  %Driving intensity
harmonicMeanG = zeros(length(sigma_vec),1);

%-----------------------
% Variable initialization
%-----------------------

w_N0 = zeros(1,N);
w_0N = zeros(1,N);
P = zeros(1,N);
A = zeros(1,N);
P_inf=1/N*ones(1,N);
for iR = 1:numRealizations
    wb=waitbar(iR/numRealizations);
    
    Energy = randn(1,N)*Delta;
    P0 = exp(-Energy*beta)/sum( exp(-Energy*beta));
    %     G0 = rand(1,N);
    G0 = linspace(0,1,N);
    randind = randperm(N);
    G0 = G0(randind);
    for iS = 1:length(sigma_vec)
        
        sigma = sigma_vec(iS);
        G = 10.^(G0*sigma - sigma/2);
        %         G = rand(1,N);
        G = G*mean(1./G);
        
        Delta_n = Energy(2:N)-Energy(1:N-1);
        Delta_n = [Delta_n, Energy(1)-Energy(N)];
        
        
        w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
        w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];
        
        w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
        w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];
        
        
        for iE = 1: length(nu_vec)
            
            w_drv = nu_vec(iE) * G * w_beta;
            w_m = (w_drv + w_bath_ccw); %w_{n-1,n}
            w_p = (w_drv + w_bath_cw);  %w_{n,n-1}
            %
            
            
            W_cl = -diag(w_drv + w_bath_cw +w_drv([N,1:N-1])+ w_bath_ccw([N,1:N-1]),0)+...
                diag(w_drv(1:N-1) + w_bath_ccw(1:N-1),1) + ...
                diag(w_drv(1:N-1) + w_bath_cw(1:N-1),-1);
            
            W_cl(N,1) = w_bath_ccw(N) + w_drv(N);
            W_cl(1,N) = w_bath_cw(N) + w_drv(N);
            
            I_mat = -diag(w_m([N,1:N-1]) - w_p,0)/N;
            
            
            B = [W_cl;ones(1,N)];
            C = zeros(N+1,1);C(N+1)=1;
            [P_cl,R] = linsolve(B,C);
            %
            
            [P_an1,P_an2,I1(iE),er1(iE)] = ness(w_m,w_p);
            
            SMF(iE) = sum(log(w_p./w_m));
            figure(1);plot(1:N,[P_cl(1),fliplr(P_cl(2:N)')],...
                1:N,(P_an1),...
                1:N,P_an2,'--');
                legend(['MEQ';'P_n';'I/G'])
        end
    end
end

