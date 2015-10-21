function [muQ,sigQ,muX,sigX,t_end] = Gillespie3(w_p,w_m)
N=length(w_p);
numRuns = 1e3;


% t_end=1/min([w_p,w_m]) *1e2;
r = w_p + w_m([N,1:N-1]);
t_end = max(1./r)*1e3;
% t_end=500;
Q=zeros(numRuns,1);
Q_prev  = zeros(numRuns,1);
H=zeros(numRuns,1);
H_prev  = zeros(numRuns,1);

X  = zeros(numRuns,1);
numTimeSteps = zeros(numRuns,1);
endTime= zeros(numRuns,1);
% site_t = zeros(numTimeSteps,numRuns);
% t_vec = zeros(numTimeSteps,numRuns);
for iR=1:numRuns
    wb=waitbar(iR/numRuns);
    
    t=0;
    site = 0;
    site_prev=0;
    iT=0;
    %     tic;
    while t<t_end
        iT=iT+1;
        n = mod(site,N)+1;
        %         n_1 = mod(site-1,N)+1;
        
        a = rand(2,1);
        
        delta_t = -1/r(n)*log(a(1));
        
        if(a(2)<w_p(n)/r(n))
            site=site+1;
        else
            site=site-1;
        end
        
        
        
        if(mod(site_prev,N)+1==3 && mod(site,N)+1==4)
            Q(iR) = Q_prev(iR) +1;
            %             H(iR) = H_prev(iR) +Delta_n(mod(site_prev,N));
        elseif(mod(site_prev,N)+1==4 && mod(site,N)+1==3)
            Q(iR) = Q_prev(iR) -1;
            %             H(iR) = H_prev(iR) -Delta_n(mod(site_prev,N));
            
        else
            Q(iR) = Q_prev(iR);
            %             H(iR) = H_prev(iR);
        end
        
        Q_prev(iR) = Q(iR);
        %         H_prev(iR) = H(iR);
        site_prev = site;
        
        t = t + delta_t;
        numTimeSteps(iR)=numTimeSteps(iR)+1;
        site_t(iT,iR)=site;
        t_vec(iT,iR)=t;
        Q_t(iT,iR)=Q(iR);
        %
    end
    X(iR)=site;
    endTime(iR) = t;
    
end

% Q_sample = Q;
% muQ = mean(Q_sample);
% sigQ = std(Q_sample);
mask=(t_vec<0.9*t_end).*(t_vec~=0);
iT_end= min(sum(mask,1));

dt=diff(t_vec(1:iT_end,:),1,1);
dt = mean(dt(:));

mask=(t_vec>(0.9*t_end-dt)).*(t_vec<0.9*t_end).*(t_vec~=0);
mask =(mask~=0);
meanQ = mean(Q_t(mask));
stdQ = std(Q_t(mask));

meanX = mean(site_t(mask));
stdX = std(site_t(mask));

% mu = 2*(w_p(1)*w_p(2)-w_m(1)*w_m(2))/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
% D = (2*(w_p(1)*w_p(2)+w_m(1)*w_m(2)) - mu^2)/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
% figure;
% hold on
% hist(site_t(mask),min(site_t(:)):max(site_t(:)));
% x=linspace(min(site_t(:)),max(site_t(:)),1000);
% plot(x,normpdf(x,meanX,stdX)*sum(mask(:)),'r','LineWidth',2)
% % plot(abs(ifft(fftshift(H2))),'g')
% 
% 
% 
% figure;
% hold on
% hist(Q_t(mask),min(Q_t(:)):max(Q_t(:)));
% QQ=linspace(min(Q_t(:)),max(Q_t(:)),1000);
% plot(QQ,normpdf(QQ,meanQ,stdQ)*sum(mask(:)),'r','LineWidth',2)


t_end=0.9*t_end;

% h=hist(site_t(mask),min(site_t(:)):max(site_t(:)));
%
%
% muQ_t=mean(Q_t,2);
% sigQ_t=std(Q_t,0,2);
% %
% % muQ = muQ_t(iT_end);
% % sigQ =sigQ_t(iT_end);
%
% muX_t=mean(site_t,2);
% sigX_t=std(site_t,0,2);

% muX = muX_t(iT_end);
% sigX =sigX_t(iT_end);

%%
%calculate mean and sig of Q within time steps of size dt
% [t_min,iT_min]=min(max(t_vec)); %find minimum end time
dt=diff(t_vec(1:iT_end,:),1,1);
dt = 10*mean(dt(:));
t=0;
ii=1;
%
mask=(t_vec>(t_end-dt)).*(t_vec<t_end).*(t_vec~=0);
mask =(mask~=0);
% meanQ = mean(Q_t(mask));
% stdQ = std(Q_t(mask));
while t<t_end
    mask=(t_vec>t).*(t_vec<t+dt).*(t_vec~=0);
    mask =(mask~=0);
    %     mask = [abs(diff(t_vec<t,1,1));zeros(1,numRuns)]~=0;
    
    meanQ(ii) = mean(Q_t(mask));
    stdQ(ii) = std(Q_t(mask));
    meanX(ii) = mean(site_t(mask));
    stdX(ii) = std(site_t(mask));
    t = t+dt;
    ii=ii+1;
end
tt = (0:ii-2)*dt+dt;

% w_pp=w_p;
% w_mm=w_m;
% for ii=1:N-1
%     
%     G_p = (1/w_p(1)+1/w_p(2) * (w_m(1)/w_p(1)))^(-1);
%     G_m = (1/w_m(2)+1/w_m(1) * (w_p(2)/w_m(2)))^(-1);
%     
%     w_p(2)=G_p;
%     w_m(2)=G_m;
%     w_p=w_p(2:end);
%     w_m=w_m(2:end);
% end
% 
% 
% 
% figure;
% axes('FontSize',24);
% hold on
% grid on
% plot(tt,stdX.^2,'or','LineWidth',4)
% 
% plot(tt,2*D*tt,'k','LineWidth',4);
% plot(tt,2*(G_p+G_m)*tt,'--b','LineWidth',4);
% xlabel('t');
% ylabel('Var(x)');
% legend(['Simulation      ';...
%         'Analytical      ';...
%         'Resistor network'],'Location','NorthWest');
%     
% axis tight;  
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/VarXvsTime.eps');

% 
muQ=meanQ(end);
sigQ=stdQ(end);
muX=meanX(end);
sigX=stdX(end);
% 
% SMF = 20;
% w_p = exp(SMF/(N*(N-1))*(0:N-1));
% w_m = exp(-SMF/(N*(N-1))*(0:N-1));
%    

% plot(tt,(G_p+G_m)*tt,'b','LineWidth',2);
% xlabel('t');
% ylabel('\sigma^2');
% legend(['var(x)_{sim}';...
% %         '2Dt         ';...
%         '(G_L+G_R)t  ']);

% figure;
% hold on
% plot(tt,stdQ.^2,'--r','LineWidth',2);
% 
% plot(tt,meanQ,'r','LineWidth',2);
% plot(tt,(G_p-G_m)*tt,'LineWidth',2);
% xlabel('t')
% ylabel('\mu');
% legend(['mean(x)_{sim}';...
%         '(G_L-G_R)t   ']);

%%

%
% figure;
% axes('FontSize',24);
% hold on;
% grid on
% plot((1:ii-1)*dt,meanQ,'r');
% plot((1:ii-1)*dt,stdQ.^2,'b');
% % plot([iT_end iT_end],[0 max(muQ_t)],'--k');
% xlabel('t[arbitrary units]');
% title(['SMF=',num2str(sum(log(w_p./w_m)))])
% legend(['\mu     ';...
%     '\sigma^2']);
% meanQ = meanQ(ii-2);
% stdQ = stdQ(ii-2);
%%
% figure;
% plot((1:ii-1)*dt,meanQ./((1:ii-1)*dt));
% figure;
% plot((1:ii-1)*dt,stdQ./((1:ii-1)*dt));
%
%
% figure;
% plot((1:ii-1)*dt,2*meanQ(1:ii-1)./stdQ(1:ii-1).^2);
%
