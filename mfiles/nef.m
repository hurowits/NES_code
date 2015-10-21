% Script for calculating current and statistics for nearest neighbor jump 
% process on a ring.
% Script produces the following figures of NEF
% - SMF vs. tau for several values of sigma.
% - I(tau;sigma) image
% - log(Abs(I)) vs log nu
% - (I,SMF) Scatter diagram

% -----------------------
% Running parameters
% -----------------------

numRealizations=1;
flagAnalytic=0;  %0 for numerical solution, 1 for analytic solution

% -----------------------
% Model parameters
% -----------------------

N=1000; %number of sites

Delta = 1;                             %Energy interval
w_beta = 1;1e-2;                       %Coupling to bath
beta = 0.1;                            %Inverse temperature
sigma_vec = 10;linspace(2,50,100);        %log-width of distribution
nu_vec = logspace(-5,5,1000);
Estrength_vec = logspace(-5,5,1000);  %Driving intensity
tau_vec=linspace(-0.2,1.2,5000);       %Scaled driving intensity
tau_vec=0.3;
harmonicMeanG = zeros(length(sigma_vec),1);

%-----------------------
% Variable initialization
%-----------------------
% if flagAnalytic ==1
%     P_analytical=zeros(numRealizations,length(Estrength_vec),N);
%     CurrentAnalytical = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
% end

w_N0 = zeros(1,N);
w_0N = zeros(1,N);
CurrentEnsemble = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
maxI = zeros(numRealizations,length(sigma_vec));
% Q_dot =  zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
P = zeros(1,N);
P_Numerical= zeros(numRealizations,length(Estrength_vec),N);
P_nu = zeros(length(nu_vec),N);
SMFEnsemble = zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
U=zeros(length(Estrength_vec),N);
Barrier = zeros(numRealizations,1);
maxV = zeros(numRealizations,1);
minV = zeros(numRealizations,1);
% rmsV = zeros(numRealizations,1);
% maxU = zeros(numRealizations,1);
% maxU1 = zeros(numRealizations,1);
% minU = zeros(numRealizations,1);
% meanU = zeros(numRealizations,1);

%
for iR = 1:numRealizations
    wb=waitbar(iR/numRealizations);
    
    Energy = randn(1,N)*Delta;
    P0 = exp(-Energy*beta)/sum( exp(-Energy*beta));
    %         G0 = rand(1,N);
    G0 = linspace(0,1,N);
%     randind = randperm(N);
    G0 = G0(randind);
    
    for iS = 1:length(sigma_vec)
        
        sigma = sigma_vec(iS);
        G = 10.^(G0*sigma - sigma/2);
        G = G*mean(1./G);
        
        nu = 1/max(G)*10.^(sigma*tau_vec);

%         engineer transition rates for increased SMF
                 for p=1:N-1
                    if (log10(G(p)*nu_vec(1))<-sigma/3)
                    
                        E1 = max(Energy(p:p+1));
                        E2 = min(Energy(p:p+1));
                        Energy(p)=E1;
                        Energy(p+1)=E2;
                    end
                end
        
        Delta_n = Energy(2:N)-Energy(1:N-1);
        Delta_n = [Delta_n, Energy(1)-Energy(N)];
        
        
        w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
        w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];
        
        w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
        w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];
        
        
        I_0 = -beta/N*dot(G,Delta_n)*nu(1) * w_beta;
        
        for iE = 1: length(tau_vec)
%             for iE = 1: length(nu_vec)
            w_drv = nu(iE) * G * w_beta;
%             w_drv = Estrength_vec(iE)^2 * G * w_beta;
            w_m = (w_drv + w_bath_ccw); %w_{n-1,n}
            w_p = (w_drv + w_bath_cw);  %w_{n,n-1}
            
            if(flagAnalytic==0)
                
                I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
                I_0 = I;
                P(1) = 1;
                for q = 1:N
                    P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
                end
                %
                er(iE)=abs(P(end)-P(1)); %error in P
                
                P = P(1:end-1)./sum(P(1:end-1)); %Normalize probabilities
                P_nu(iE,:) = P;


                %                 Q_dot(iR,iE,iS) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
                %                     - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ; %equation (31)
                %                 CurrentEnsemble(iR,iE,iS) = I;
                                P_Numerical(iR,iE,:) = P;
                CurrentEnsemble(iR,iE,iS) =(P(1)*w_p(1)-P(2)*w_m(1))*N ; %adjust current according to normalization
                
                %                 % Find first maximum of current vs. nu
                %                 if (abs(CurrentEnsemble(iR,iE,iS))> abs(CurrentEnsemble(iR,iE-1,iS)))
                %                     maxI(iR,iS) = CurrentEnsemble(iR,iE,iS);
                %                 else
                %                     break
                %                 end
                %
                SMFEnsemble(iR,iE,iS) = sum(log(w_m)-log(w_p));
                V1 = [0 cumsum(log(w_m./w_p))];
                Vflat = V1-V1(end)*(0:N)/(N);
                Vflat =  Vflat(1:end-1);
                Barrier(iR,iE) = max(Vflat)-min(Vflat);
                maxV(iR,iE) = max(Vflat);
                                U(iE,:)=Vflat;
                                V(iE,:)=V1;
                
                %                 minV(iR,iE) = min(Vflat);
                %                 rmsV(iR,iE) = std(Vflat);
                %                 meanU(iR,iE)= mean(Vflat);
                %                 minU(iR,iE) =  min(Vflat - mean(Vflat));
                %                 maxU(iR,iE)= max(Vflat - mean(Vflat));
                %                 maxU1(iR,iE)= max(Vflat)- mean(Vflat);
                
            elseif(flagAnalytic==1)
                
                for l = 1:N
                    
                    V1 = [0 cumsum(log(w_m./w_p))];                 %V1(n)=V(n->0), V1(end) = -SMF
                    V2 = [fliplr(cumsum(log(fliplr(w_p./w_m)))) 0]; %V2(n)=V(n->N), V2(1) = SMF
                    
                    w_N0(l) = 1./sum(exp(V1(1:end-1))./w_p);        %exp(V(n-1->0))/w_{n,n-1}
                    w_0N(l) = 1./sum(exp(V2(2:end))./w_m);          %exp(V(n->N))/w_{n-1,n}
                    
                    w_p = circshift(w_p,[0,1]);
                    w_m = circshift(w_m,[0,1]);
                    
                end
                
                CurrentAnalytical(iR,iE) = 1./sum(1./(w_N0-w_0N));
                P_analytical(iR,iE,:) = CurrentAnalytical(iR,iE)./(w_N0-w_0N);
            end
            
        end
    end
end
max(er)
% B = diff(log(P_nu))./repmat(diff(nu_vec'),[1 N]);
chi = (B(:,1)*P0(1) - B(:,2)*P0(2))*G(1);
chi2 = (B(:,2)*P0(2) - B(:,3)*P0(3))*G(2);
%%

