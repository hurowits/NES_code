%%
clear all
N=1e2;
numRealizations=1;
% numRealizations=1e4;
randInd1 = randperm(N);
randInd2 = randperm(N);
systemModel = 1; %0 for chain, 1 for ring
flagAnalytic=0;  %0 for numerical solution, 1 for analytic solution


Delta = 1;
w_beta = 1;
beta = 0.1;
sigma_vec = 10;[1,2,4];[1,3,5,7];1:10;
sigma = sigma_vec(1);
% Estrength_vec = logspace(-5,3,100);
Estrength_vec = logspace(-5,3,1000);
% Estrength_vec = [1e-2, 1, 1e3]
% Estrength_vec = 500;
Current  = zeros(numRealizations, length(Estrength_vec));
CurrentAnalytical =  zeros(numRealizations, length(Estrength_vec));
Q_dot =  zeros(numRealizations, length(Estrength_vec));
SMF =  zeros(numRealizations, length(Estrength_vec));
P = zeros(1,N);
% Estrength_vec = [logspace(-3,-2,10), linspace(10^-2,10^2,80), logspace(2,3,10)];
% Estrength_vec = linspace(0,1e6,100);

SMF = zeros(length(Estrength_vec), length(sigma_vec));
seed=1;
% G =genRand2(sigma_vec,randperm(N),N,'sigma');
%
% G = sigma_vec*rand(1,N);
% G = exp(G);
G = logspace(-sigma_vec,sigma_vec,N);
% G = genRand(sigma,randperm(N),N+1,'sigma');


% rndNumNrg = rand(1,N);
% Energy = rand(1,N)*Delta - Delta/2;
Energy = linspace(-Delta/2,Delta/2,N);
Energy = Energy(randInd2);
Delta_n = Energy(2:N)-Energy(1:N-1);
Delta_n = [Delta_n, Energy(1)-Energy(N)];


w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];

w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];


for iR = 1:numRealizations
%     G = G(randperm(N));
    G = 10.^(rand(1,N)*sigma - sigma/2);
    G = G*mean(1./G);
    
    
    
    
    I_0=0;
    I_0 = -beta/N*dot(G,Delta_n)*Estrength_vec(1)^2;
    %     I_0=1e-10;
    % I_0=0;
    
    for iE = 1: length(Estrength_vec)
            wb=waitbar(iE/length(Estrength_vec));

        w_drv = Estrength_vec(iE)^2 * G * w_beta;
        w_m = (w_drv + w_bath_ccw);%/max(w_drv);
        w_p = (w_drv + w_bath_cw);%/max(w_drv);
        
        if(systemModel==1)
%             SMF(iE,iSigma) = sum(log(w_bath_ccw+w_drv)-log(w_bath_cw+w_drv))/beta;
            SMF(iR,iE) = sum(log(w_m)-log(w_p))/beta;
            I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
            I_0 = I;
            P(1) = 1;
            for q = 1:N
                P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
            end
            
            P = P(1:end-1)./sum(P(1:end-1));
            Q_dot(iR,iE,systemModel) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ; %equation (31)
            Current(iR,iE) = I;
%             

            if(flagAnalytic==1)
            den=0;
            ;
                
            for t = 1:N
                for s=1:N
                    den = den + prod(w_p(s+1:N))*prod(w_m(s:-1:2));
                end
                w_p = circshift(w_p,[0,1]);
                w_m = circshift(w_m,[0,1]);
            end
            
            A(iR) = prod(w_p);
            B(iR) = prod(w_m);
            D(iR) = den;
            CurrentAnalytical(iE) = N*(prod(w_p) - prod(w_m))/den;
            end
%             %equation (25)
            %                 landscape = cumsum((log(w_m)-log(w_p))/beta,2);
            %                 maxVar = max(landscape)-min(landscape);
            %                 Current_Estimate(iE,iSigma) = sum(1./(w_m+w_p))^-1* exp(-beta*maxVar/2)*sinh(beta*SMF(iE,iSigma)/2);
            
            %         Current_Estimate2(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-2*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         Current_Estimate3(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-3*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         Current_Estimate4(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-4*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         w_harm(iE,iSigma) = mean(2./(w_m+w_p));
        end
    end
   
end

close(wb);
%%
% sigma = std(Current,1);
% mu = mean(Current,1);
% [h,x]=hist(Current,100);
% figure;plot(x,h/numRealizations);
% figure;hold on
% plot(x,normpdf(x,mu(1),sigma(1)));
% plot(x,normpdf(x,mu(2),sigma(2)));
% plot(x,normpdf(x,mu(3),sigma(3)));

% plot(x,1./sqrt(2*pi*sigma(2).^2).*exp(-(x-mu(2)).^2./(2*sigma(2).^2)));
% plot(x,1./sqrt(2*pi*sigma(3).^2).*exp(-(x-mu(3)).^2./(2*sigma(3).^2)));

% xlabel('I');
% ylabel('P(I)');
% legend(num2str((Estrength_vec').^2))
%% SMF and Current log-log
 %% SMF and Current log-log
    iSigma=1;
    figure(1);
    %hold on
    %axes('FontSize',14);
    loglog(Estrength_vec.^2,abs(SMF(1,:)),'-r', 'LineWidth',2);hold on
    loglog(Estrength_vec.^2,abs(Current(1,:)),'b', 'LineWidth',2);
    loglog(Estrength_vec.^2,abs(CurrentAnalytical(1,:)),':r', 'LineWidth',2);

%     loglog(Estrength_vec.^2,abs(Current_Estimate),'b--', 'LineWidth',2);hold on
    % loglog(Estrength_vec.^2,abs(1./w_harm(:,iSigma)),'--r', 'LineWidth',2);
    % loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_0(:,iSigma)),'--r', 'LineWidth',2);
    % loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_inf(:,iSigma)),'--r', 'LineWidth',2);
    SinaiLimit = sqrt(N)*Delta/2;
    loglog(Estrength_vec.^2, sqrt(N)*Delta/2.*ones(length(Estrength_vec),1),'-k', 'LineWidth',2);
    % I_typ = std(G)/sqrt(N) *Delta*beta*Estrength_vec.^2;
    I_0 = -beta/N*dot(G(iSigma,:),Delta_n)*Estrength_vec.^2;
    I_inf = beta/N*dot(1./G(iSigma,:),Delta_n)*ones(length(Estrength_vec),1)/mean(1./G(iSigma,:));
    loglog(Estrength_vec.^2,abs(I_0),'--g', 'LineWidth',2);
    loglog(Estrength_vec.^2,abs(I_inf),'--g','LineWidth',2);
    loglog(Estrength_vec.^2,2/beta*ones(size(Estrength_vec)),'--m','LineWidth',2);
    % loglog(Estrength_vec.^2,abs(I_typ),'--g');
    % axis([Estrength_vec(1).^2 Estrength_vec(end).^2
    % min(abs(Current(:,iSigma))), 10*SinaiLimit]);
    % axis([1e-5 1e5 1e-5 1e1]);
    currentLims = linspace(min(abs(Current(1,:))),max(abs(SMF(1,:))));
%     axis([1e-10 Est rength_vec(end).^2 currentLims(1), 10*SinaiLimit]);
    
    loglog(ones(length(currentLims),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
    loglog(ones(length(currentLims),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
%     loglog(ones(length(currentLims),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
%     loglog(ones(length(currentLims),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
    
    legend(['|SMF|    ';...
        '|Current|']);
    xlabel('\epsilon^2');
    ylabel('|Current| and |SMF|');
    grid on
    hold off