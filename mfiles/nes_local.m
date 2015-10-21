% clear all
N=1000; %number of sites
numRealizations=10000;
% numRealizations=1
flagAnalytic=1;  %0 for numerical solution, 1 for analytic solution


Delta = 1;
w_beta = 1e-2;
beta = 0.1;
% sigma_vec =linspace(0,10,50);
% Estrength_vec = logspace(-5,3,1000);

% Estrength_vec = Estrength_vec(601);
% sigma_vec = sigma_vec([30]);
sigma_vec=6;
Estrength_vec=1;
harmonicMeanG = zeros(length(sigma_vec),1);
% Current  = zeros( length(Estrength_vec),length(sigma_vec));
% Current_Estimate  = zeros( length(Estrength_vec),length(sigma_vec));
% CurrentAnalytical =  zeros( length(Estrength_vec),length(sigma_vec));
CurrentAnalytical =  zeros( numRealizations);
A =  zeros( numRealizations,1);
B =  zeros( numRealizations,1);
D =  zeros( numRealizations,1);
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
Energy = randn(1,N)*Delta;
Delta_n = Energy(2:N)-Energy(1:N-1);
Delta_n = [Delta_n, Energy(1)-Energy(N)];


w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];

w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];


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
            w_m = (w_drv + w_bath_ccw);
            w_p = (w_drv + w_bath_cw);
            %             if(flagAnalytic==0)
            SMFEnsemble(iR,iE,iS) = sum(log(w_m)-log(w_p));
            I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
            I_0 = I;
            P(1) = 1;
            for q = 1:N
                P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
            end
            
            P = P(1:end-1)./sum(P(1:end-1));
            
            Q_dot(iR,iE,iS) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ; %equation (31)
            CurrentEnsemble(iR,iE,iS) = I;
            %
            if(flagAnalytic==1)
%                 den=0;
%                 for t = 1:N
%                     for s=1:N
%                         den = den + prod(w_p(s+1:N))*prod(w_m(s:-1:2));
%                     end
%                     w_p = circshift(w_p,[0,1]);
%                     w_m = circshift(w_m,[0,1]);
%                 end
                
%                 A(iR) = prod(w_p);
%                 B(iR) = prod(w_m);

                A(iR) = sum(log(w_p./w_m));
                B(iR) = sum(log(w_p.*w_m));
                
               
%                 D(iR) = den;
                
%                 A(iR) = exp( sum(log(w_p)) - log(den));
%                 B(iR) = exp( sum(log(w_m)) - log(den));
%                 CurrentAnalytical(iR) = A(iR)-B(iR);%N*(prod(w_p) - prod(w_m))/den;
%                 w_drv2=w_drv;
%                 for p=1:N
%                     temp = w_drv2(p);
%                     w_drv2(p)=Delta_n(p)*beta/2;
%                     numeratorEst(iR) = numeratorEst(iR) + prod(w_drv2);
%                     w_drv2(p)=temp;
%                 end
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
 X = sqrt(exp(A).*exp(B));
 Y = sqrt(exp(B)./exp(A));
1
%% SMF and Current log-log
iSigma=1;
figure(1);
%hold on
%axes('FontSize',14);
loglog(Estrength_vec.^2,abs(SMF(:,1)),'-r', 'LineWidth',2);hold on
loglog(Estrength_vec.^2,abs(Current(:,1)),'b', 'LineWidth',2);
% loglog(Estrength_vec.^2,abs(CurrentAnalytical(1,:,2)),':r', 'LineWidth',2);

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
%
% %% SMF and Current log-log
% figure(2);
% %hold on
% %axes('FontSize',14);
% semilogy(sigma_vec,abs(SMF(1,:)),'-r', 'LineWidth',2);hold on
% semilogy(sigma_vec,abs(Current(1,:)),'b', 'LineWidth',2);
% % loglog(Estrength_vec.^2,abs(CurrentAnalytical(1,:,2)),':r', 'LineWidth',2);
%
% %     loglog(Estrength_vec.^2,abs(Current_Estimate),'b--', 'LineWidth',2);hold on
% % loglog(Estrength_vec.^2,abs(1./w_harm(:,iSigma)),'--r', 'LineWidth',2);
% % loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_0(:,iSigma)),'--r', 'LineWidth',2);
% % loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_inf(:,iSigma)),'--r', 'LineWidth',2);
% % SinaiLimit = sqrt(N)*Delta/2;
% % loglog(Estrength_vec.^2, sqrt(N)*Delta/2.*ones(length(Estrength_vec),1),'-k', 'LineWidth',2);
% % % I_typ = std(G)/sqrt(N) *Delta*beta*Estrength_vec.^2;
% % I_0 = -beta/N*dot(G(iSigma,:),Delta_n)*Estrength_vec.^2;
% % I_inf = beta/N*dot(1./G(iSigma,:),Delta_n)*ones(length(Estrength_vec),1)/mean(1./G(iSigma,:));
% % loglog(Estrength_vec.^2,abs(I_0),'--g', 'LineWidth',2);
% % loglog(Estrength_vec.^2,abs(I_inf),'--g','LineWidth',2);
% % loglog(Estrength_vec.^2,2/beta*ones(size(Estrength_vec)),'--m','LineWidth',2);
% % % loglog(Estrength_vec.^2,abs(I_typ),'--g');
% % % axis([Estrength_vec(1).^2 Estrength_vec(end).^2
% % % min(abs(Current(:,iSigma))), 10*SinaiLimit]);
% % % axis([1e-5 1e5 1e-5 1e1]);
% % currentLims = linspace(min(abs(Current(1,:))),max(abs(SMF(1,:))));
% % %     axis([1e-10 Est rength_vec(end).^2 currentLims(1), 10*SinaiLimit]);
% %
% % loglog(ones(length(currentLims),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
% % loglog(ones(length(currentLims),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
% %     loglog(ones(length(currentLims),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
% %     loglog(ones(length(currentLims),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
%
% legend(['|SMF|    ';...
%     '|Current|']);
% xlabel('\sigma');
% ylabel('|Current| and |SMF|');
% grid on
% hold off
