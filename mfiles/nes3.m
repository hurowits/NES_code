% clear all
N=10000; %number of sites
% numRealizations=1000;
numRealizations=1;
flagAnalytic=1;  %0 for numerical solution, 1 for analytic solution

mode = 2;% 1=chain, 2=ring

Delta = 1;
w_beta = 1;1e-2;
beta = 0.1;
% sigma_vec =linspace(0,10,50);
Estrength_vec = logspace(-5,3,100);

% Estrength_vec = Estrength_vec(601);
% sigma_vec = sigma_vec([30]);
sigma_vec=6;
Estrength_vec=1;
harmonicMeanG = zeros(length(sigma_vec),1);
% Current  = zeros( length(Estrength_vec),length(sigma_vec));
% Current_Estimate  = zeros( length(Estrength_vec),length(sigma_vec));
% CurrentAnalytical =  zeros( length(Estrength_vec),length(sigma_vec));
P_analytical=zeros(numRealizations,length(Estrength_vec),N);
CurrentAnalytical = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
w_N0 = zeros(1,N);
w_0N = zeros(1,N);
% numeratorEst = zeros(numRealizations,1);
CurrentEnsemble = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
% CurrentEstimateEnsemble = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
% wAverageEnsemble = zeros(numRealizations, length(Estrength_vec),length(sigma_vec));
Q_dot =  zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
P = zeros(1,N);

% SMF = zeros(length(Estrength_vec), length(sigma_vec));
SMFEnsemble = zeros(numRealizations,length(Estrength_vec),length(sigma_vec));

% landscape = zeros(length(Estrength_vec), length(sigma_vec),N);
% landscape2 = zeros(numRealizations,length(Estrength_vec), length(sigma_vec),N);
% barrier = zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
% maxima = zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
% minima = zeros(numRealizations,length(Estrength_vec),length(sigma_vec));
% seed=1;

% Energy = linspace(-Delta/2,Delta/2,N);
% Energy = Energy(randInd2);
if (mode==1)
    Energy = randn(1,N+1)*Delta;
    Delta_n = Energy(2:N+1)-Energy(1:N);
    
    
    w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N)-Energy(2:N+1)))); %w_plus
    
    
    w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N)-Energy(2:N+1)))); %w_minus
    
else
    Energy = randn(1,N)*Delta;
    Delta_n = Energy(2:N)-Energy(1:N-1);
    Delta_n = [Delta_n, Energy(1)-Energy(N)];
    
    
    w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
    w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];
    
    w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
    w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];
end
for iR = 1:numRealizations
    wb=waitbar(iR/numRealizations);
    %     tic;
    for iS = 1:length(sigma_vec)
        sigma = sigma_vec(iS);
        G = 10.^(rand(1,N)*sigma - sigma/2);
        %         harmonicMeanG(iS) = mean(1./G)^-1;
        G = G*mean(1./G);
        
        %     I_0=0;
        I_0 = -beta/N*dot(G,Delta_n)*Estrength_vec(1)^2 * w_beta;
        
        for iE = 1: length(Estrength_vec)
            
            w_drv = Estrength_vec(iE)^2 * G * w_beta;
            w_m = (w_drv + w_bath_ccw); %w_{n,n-1}
            w_p = (w_drv + w_bath_cw);  %w_{n-1,n}
            %             if(flagAnalytic==0)
            SMFEnsemble(iR,iE,iS) = sum(log(w_m)-log(w_p));
            I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
            I_0 = I;
            P(1) = 1;
            for q = 1:N
                P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
            end
            
            P = P(1:end-1)./sum(P(1:end-1));
            %
            %             Q_dot(iR,iE,iS) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
            %                 - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ; %equation (31)
            CurrentEnsemble(iR,iE,iS) = I;
            
            if(flagAnalytic==1)
                
                for l = 1:N
                    
                    V1 = [0 cumsum(log(w_p./w_m))];                 %V1(n)=V(n->0), V1(end) = -SMF
                    V2 = [fliplr(cumsum(log(fliplr(w_m./w_p)))) 0]; %V2(n)=V(n->N), V2(1) = SMF
                    
                    w_N0(l) = 1./sum(exp(V1(1:end-1))./w_m);        %exp(V(n-1->0))/w_{n,n-1}
                    w_0N(l) = 1./sum(exp(V2(2:end))./w_p);          %exp(V(n->N))/w_{n-1,n}
                    %
                    w_p = circshift(w_p,[0,1]);
                    w_m = circshift(w_m,[0,1]);
                    
                end
                CurrentAnalytical(iR,iE) = 1./sum(1./(w_N0-w_0N));
                P_analytical(iR,iE,:) = CurrentAnalytical(iR,iE)./(w_N0-w_0N);
            end
            %equation (25)
            %             landscape(iE,iS,:) = cumsum((log(w_m)-log(w_p))/beta,2);
            %             landscape2(iR,iE,iS,:) = cumsum((log(w_m)-log(w_p))/beta,2);
            % maxVar = max(landscape(iE,iS,:))-min(landscape(iE,iS,:));
            % barrier(iR,iE,iS) = maxVar;
            % maxima(iR,iE,iS) = max(landscape(iE,iS,:));
            % minima(iR,iE,iS) = min(landscape(iE,iS,:));
            % CurrentEstimateEnsemble(iR,iE,iS,:) = 2*sum(1./(w_m+w_p))^-1* exp(-beta*maxVar/2)*sinh(beta*SMF(iE,iS)/2);
            % SMFEnsemble(iR,iE,iS)=SMF;
            %             wAverageEnsemble(iR,iE,iS)=sum(1./(w_m+w_p))^-1;
            % Current_Estimate2(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-2*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            % Current_Estimate3(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-3*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            % Current_Estimate4(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-4*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            % w_harm(iE,iSigma) = mean(2./(w_m+w_p));
            
        end
    end
    
    %     save (['datasets/prcs',num2str(processNumber),'/nesDataset',num2str(iR),'.mat']);
    %     CurrentEnsemble(iR,:,:,:) = Current;
    
    %     close(wb);
    %     toc;
end
%%
% figure;
% plot(log10(Estrength_vec.^2),(abs(CurrentAnalytical*N)),'b');
% hold on
% plot(log10(Estrength_vec.^2),(abs(CurrentEnsemble)),'r');
% legend(['analytical';...
%     'numerics  ']);
% 
% figure;
% plot(sort(abs(CurrentAnalytical*N)),log10(numRealizations:-1:1),'b');
% hold on
% plot(sort(abs(CurrentEnsemble)),log10(numRealizations:-1:1),'r');

figure;
axes('FontSize',24);
hold on;
grid on
plot(1:N, V1(1:end-1),'b');
plot(1:N,squeeze(P_analytical(1,1,end:-1:1))*N,'g');
plot(1:N,N*exp(-V1(2:end))./sum(exp(-V1(2:end))),'r');
legend(['V(x)      ';...
    'N*P(x)    ';...
    'exp(-V(x))'])
xlabel('bond');