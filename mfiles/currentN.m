%%
clear all
N=1e5;
randInd1 = randperm(N);
randInd2 = randperm(N);
systemModel = 1; %0 for chain, 1 for ring


Delta = 1;
w_beta = 1;
beta = 0.1;
sigma_vec = 4;[1,2,4];[1,3,5,7];1:10;
Estrength_vec = logspace(-2,3,100);
Estrength_vec = logspace(-5,3,1000);
Q_dot = zeros(size(Estrength_vec));
P = zeros(1,N);
landscape = zeros(length(Estrength_vec),N);
% Estrength_vec = [logspace(-3,-2,10), linspace(10^-2,10^2,80), logspace(2,3,10)];
% Estrength_vec = linspace(0,1e6,100);

SMF = zeros(length(Estrength_vec), length(sigma_vec));
seed=1;
% G =genRand2(sigma_vec,randperm(N),N,'sigma');
% 
% G = sigma_vec*rand(1,N);
% G = exp(G);
G = logspace(-sigma_vec,sigma_vec,N);
G = G(randInd1);
G = G*mean(1./G);

%
% G = logspace(-3,0,N);
% G = G*mean(1./G);
% G = G(randperm(N));

w_drv_0 = G*w_beta ;

% rndNumNrg = rand(1,N);
% Energy = rand(1,N)*Delta - Delta/2;
Energy = linspace(-Delta/2,Delta/2,N);
Energy = Energy(randInd2);
Delta_n = Energy(2:N)-Energy(1:N-1);
Delta_n = [Delta_n, Energy(1)-Energy(N)];

% Energy = [1,3,5,10];

Delta_n = -[Energy(1:N-1) - Energy(2:N) Energy(N)-Energy(1)];

% V_0 =[0.1,3,10];
%V_23 V_13 V_12
% G = [10,0.1,3,15].^2;

w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];

w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];




for systemModel =1
    
    
    for iSigma = 1:length(sigma_vec)
        I_0=0;
        I_0 = -beta/N*dot(G(iSigma,:),Delta_n)*Estrength_vec(1)^2;
        %     I_0=1e-10;
        % I_0=0;
        for iE = 1: length(Estrength_vec)
            wb=waitbar(iE/length(Estrength_vec));
            w_drv = Estrength_vec(iE)^2 * G(iSigma,:) * w_beta;
            w_m = w_drv + w_bath_ccw;
            w_p = w_drv + w_bath_cw;
            if(systemModel==1)
                %         SMF(iE,iSigma) =
                %         sum(log(w_bath_ccw+w_drv)-log(w_bath_cw+w_drv))/beta;
                SMF(iE,iSigma) = sum(log(w_m)-log(w_p))/beta;
%                 I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
%                 I_0 = I;
%                 P(1) = 1;
%                 for q = 1:N
%                     P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
%                 end
%                 
%                 P = P(1:end-1)./sum(P(1:end-1));
%                 Q_dot(iE,systemModel) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
%                     - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ; %equation (31)
%                 Current(iE,iSigma) = I;
                %equation (25)
                landscape(iE,:) = cumsum((log(w_m)-log(w_p))/beta,2);
%                 maxVar = max(landscape(iE,:))-min(landscape(iE,:));
%                 Current_Estimate(iE,iSigma) = sum(1./(w_m+w_p))^-1* exp(-beta*maxVar/2)*sinh(beta*SMF(iE,iSigma)/2);
            
            else
                
                w_m(N)=0;
                w_p(N)=0;
                
                P(1) = 1;
                for q = 1:N
                    P(q+1) = (w_p(q)*P(q)) / w_m(q);
                end
                
                P = P(1:end-1)./sum(P(1:end-1));
                Q_dot(iE,systemModel) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P'...
                    - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P([N,1:N-1])' ;
            end
            %         Current_Estimate2(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-2*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         Current_Estimate3(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-3*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         Current_Estimate4(iE,iSigma) = 1/N*mean(1./(w_m+w_p))^-1* exp(-4*beta*SMF(iE,iSigma)/2)*sinh(beta*SMF(iE,iSigma)/2);
            %         w_harm(iE,iSigma) = mean(2./(w_m+w_p));
        end
    end
end
close(wb);


%% SMF and Current log-log
iSigma=1;
figure;
%hold on
%axes('FontSize',14);
loglog(Estrength_vec.^2,abs(SMF(:,iSigma)),'-r', 'LineWidth',2);hold on
loglog(Estrength_vec.^2,abs(Current(:,iSigma)),'b', 'LineWidth',2);
loglog(Estrength_vec.^2,abs(Current_Estimate),'b--', 'LineWidth',2);hold on
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
currentLims = linspace(min(abs(Current)),max(abs(SMF)),length(Estrength_vec));
axis([1e-10 Estrength_vec(end).^2 currentLims(1), 10*SinaiLimit]);

loglog(ones(length(Estrength_vec),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);

legend(['|SMF|    ';...
    '|Current|']);
xlabel('\epsilon^2');
ylabel('|Current| and |SMF|');
grid on
%% Current semi-logx
iSigma=1;
figure;
axes('FontSize',14);
% semilogx(Estrength_vec.^2,abs(SMF(:,iSigma)),'-r', 'LineWidth',2);hold on
loglog(Estrength_vec.^2,abs(Current(:,iSigma)),'b', 'LineWidth',2);hold on
% loglog(Estrength_vec.^2,abs(Current_Estimate),'b--', 'LineWidth',2);hold on

SinaiLimit = sqrt(N)*Delta/2;

% I_typ = std(G)/sqrt(N) *Delta*beta*Estrength_vec.^2;
I_0 = -beta/N*dot(G(iSigma,:),Delta_n)*Estrength_vec.^2;
I_inf = beta/N*dot(1./G(iSigma,:),Delta_n)*ones(length(Estrength_vec),1)/mean(1./G(iSigma,:));
% semilogx(Estrength_vec.^2,abs(I_0),'--g', 'LineWidth',2);
semilogx(Estrength_vec.^2,abs(I_inf),'--g','LineWidth',2);
% semilogx(Estrength_vec.^2,2/beta*ones(size(Estrength_vec)),'--g','LineWidth',2);
% loglog(Estrength_vec.^2,abs(I_typ),'--g');
% axis([Estrength_vec(1).^2 Estrength_vec(end).^2
% min(abs(Current(:,iSigma))), 10*SinaiLimit]);
% axis([1e-10 1e8 0 5*1e-4]);

% axis([1e-6 Estrength_vec(end).^2 currentLims(1), currentLims(end)]);
% %
% semilogx(ones(length(Estrength_vec),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
% semilogx(ones(length(Estrength_vec),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
% % semilogx(ones(length(Estrength_vec),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
% semilogx(ones(length(Estrength_vec),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);

legend(['|Current|']);
xlabel('\epsilon^2');
ylabel('|Current|');
grid on
%% SMF semi-logx
iSigma=1;
figure;
axes('FontSize',14);
semilogx(Estrength_vec.^2,(SMF(:,iSigma)),'-r', 'LineWidth',2);hold on
semilogx(Estrength_vec.^2, sqrt(N)*Delta/2.*ones(length(Estrength_vec),1),'-k', 'LineWidth',2);
% axis([Estrength_vec(1).^2 Estrength_vec(end).^2
% min(abs(Current(:,iSigma))), 10*SinaiLimit]);
% axis([1e-10 1e8 0 5*1e-4]);
semilogx(ones(length(Estrength_vec),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
semilogx(ones(length(Estrength_vec),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',2);
% semilogx(ones(length(Estrength_vec),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
% semilogx(ones(length(Estrength_vec),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);

legend(['SMF']);
xlabel('\epsilon^2');
ylabel('SMF');
grid on

%% current vs. smf
figure;
axes('FontSize',14);
scatter(log10(abs(SMF(1:end))),log10(abs(Current(1:end))),3, log10(Estrength_vec.^2))
xlabel('abs(SMF)');
ylabel('abs(Current)');

%% Heating rate: ring vs. chain
figure;
axes('FontSize',14);
% semilogx(Estrength_vec.^2,(Q_dot) ,'r','LineWidth',2);hold on
semilogx(Estrength_vec.^2,(Q_dot(:,1)-Q_dot(:,2)) ,'--r','LineWidth',2);hold on
% semilogx(Estrength_vec.^2,(Q_dot(:,1)) ,'--b','LineWidth',2);hold on
% semilogx(Estrength_vec.^2,(Q_dot(:,2)) ,'b','LineWidth',2);hold on
semilogx(Estrength_vec.^2,Current.*SMF,'LineWidth',2);
% semilogx(Estrength_vec.^2,-(Q_dot2'+Current.*SMF),'g','LineWidth',2);
QdotLims = linspace(min(Q_dot(:,1)-Q_dot(:,2)),max(Q_dot(:,1)-Q_dot(:,2)));
% QdotLims = linspace(min(abs(-                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Current.*SMF)),max(abs(-Current.*SMF)));
axis([Estrength_vec(1).^2 Estrength_vec(end).^2 QdotLims(1), QdotLims(end)]);
% %
loglog(ones(length(Estrength_vec),1).*1./min(G(iSigma,:),[],2),QdotLims,'k--', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*1./max(G(iSigma,:),[],2),QdotLims,'k--', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*1./mean(G(iSigma,:),2),QdotLims,'k:', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*mean(1./G(iSigma,:),2),QdotLims,'k:', 'LineWidth',2);
grid on
xlabel('\epsilon^2');
ylabel('Qdot[Ring] - Qdot[Chain]');
legend(['Qdot[Ring] - Qdot[Chain]';...
        'Entropy Production      ']);
%% Heating rate and entropy production
figure;
Q_dot_inf = beta/N * sum(Delta_n.^2);
axes('FontSize',14);
loglog(Estrength_vec.^2,Q_dot(:,1) ,'r','LineWidth',2);hold on
% semilogx(Estrength_vec.^2,-Q_dot(:,2) ,'g','LineWidth',2);hold on
semilogx(Estrength_vec.^2,abs(Q_dot(:,1)-Q_dot(:,2)) ,'--r','LineWidth',2);hold on
semilogx(Estrength_vec.^2,-Current.*SMF,'LineWidth',2);
semilogx(Estrength_vec.^2,Q_dot_inf.*ones(size(Estrength_vec)) ,'--g','LineWidth',2);hold on
QdotLims = linspace(min(Q_dot(:,1)),max(Q_dot(:,1)));
loglog(ones(length(Estrength_vec),1).*1./min(G(iSigma,:),[],2),QdotLims,'k--', 'LineWidth',2);
loglog(ones(length(Estrength_vec),1).*1./max(G(iSigma,:),[],2),QdotLims,'k--', 'LineWidth',2);
% loglog(ones(length(Estrength_vec),1).*1./mean(G(iSigma,:),2),QdotLims,'k:', 'LineWidth',2);
% loglog(ones(length(Estrength_vec),1).*mean(1./G(iSigma,:),2),QdotLims,'k:', 'LineWidth',2);
grid on
xlabel('\epsilon^2','FontSize',27);
ylabel('EAR','FontSize',27);
legend(['Qdot              ';...
        'Entropy procuction']);
%% Energy landscape
a = jet;
l1 = min(landscape(:))-max(landscape(:));
l2 = max(landscape(:)) - min(landscape(:));
landscape(1) = l1;landscape(2) = l2;
figure;imagesc(landscape);colormap(a([33:64,1:32],:))