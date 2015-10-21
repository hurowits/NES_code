% function P = Gillespie(w_p,w_m)
clear site_t t_vec Q
N=length(w_p);
numRuns = 1e3;
numTimeSteps = 1e5;

t_end=1/min([w_p,w_m]);
r = w_p + w_m([N,1:N-1]);
Q=zeros(numTimeSteps,numRuns);

site_t = zeros(numTimeSteps,numRuns);
t_vec = zeros(numTimeSteps,numRuns);
for iR=1:numRuns
    wb=waitbar(iR/numRuns);
    
    t=0;
    site = 0;
    %     tic;
    for iT=2:numTimeSteps
        n = mod(site,N)+1;
        n_1 = mod(site-1,N)+1;
        
        a = rand(2,1);
        
        delta_t = -1/r(n)*log(a(1)/r(n));
        
        if(a(2)<w_p(n)/r(n))
            site=site+1;
        else
            site=site-1;
        end
        
        site_t(iT,iR) = site;%mod(site,N)+1;
        
        if(mod(site_t(iT-1,iR),N)+1==3 && mod(site_t(iT,iR),N)+1==4)
            Q(iT,iR) = Q(iT-1,iR) +1;
        elseif(mod(site_t(iT-1,iR),N)+1==4 && mod(site_t(iT,iR),N)+1==3)
            Q(iT,iR) = Q(iT-1,iR) -1;
        else
            Q(iT,iR) = Q(iT-1,iR);
        end
        
        
        t = t + delta_t;
        t_vec(iT,iR) = t;
        
    end
    
end
%%

[t_min,iT_min]=min(max(t_vec)); %find minimum end time

inds = sum(t_vec<t_min,1); %index of minimum end time of each realization
figure;
axes('FontSize',24);

plot(t_vec(t_vec(:,1)<t_min,1),site_t(t_vec(:,1)<t_min,1),'b');
hold on
plot(t_vec(t_vec(:,1)<t_min,1),Q(t_vec(:,1)<t_min,1)*N,'r');
% plot(t_vec(t_vec(:,1)<t_min,1),N* ceil(site_t(t_vec(:,1)<t_min,1)/N),'--g');
xlabel('t');
ylabel('x(t),Q(t)');
axis ([0 t_min min(site_t(t_vec(:,1)<t_min)) 1.1*max(site_t(t_vec(:,1)<t_min))])
legend(['x(t)';'Q(t)'])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/x_t_Sat.eps');

%%
mask = [abs(diff(t_vec<t_min,1,1));zeros(1,numRuns)]~=0;

inds = sum(t_vec<t_min);
h=hist(mod(site_t(mask),N)+1,1:N);

figure;
axes('FontSize',22);
stem(1:N,h, 'LineWidth',4)
xlabel('n');
ylabel('P(n)');
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/p_n_Sat.eps');

figure;
axes('FontSize',24);

Q_sample = Q(mask);
mu1 = mean(Q_sample);
sig1 = std(Q_sample);
[h,x] = hist(Q_sample, -max(abs(Q_sample)):1:max(abs(Q_sample)));
%  hist(Q_sample, -max(abs(Q_sample)):1:max(abs(Q_sample)));
hold on
plot(x,normpdf(x,mu1,sig1)*length(Q_sample),'r','LineWidth',4)
stem(x,h,'MarkerSize',5)
xlabel('Q');
ylabel('P(Q)');
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/p_Q_Sat.eps');

% %%
% dt=diff(t_vec,1,1);
% dt = 100*mean(dt(:));
% ii=0;
% t=0;
% meanQ=zeros(1,length(0:dt:t_min));
% stdQ=zeros(1,length(0:dt:t_min));
% t_vec2=zeros(1,length(0:dt:t_min));
% while t<t_min
%     ii=ii+1    ;
%     t = t+dt;
%     t_vec2(ii)=t;
%     mask = [abs(diff(t_vec<t,1,1));zeros(1,numRuns)]~=0;
%     meanQ(ii) = mean(Q(mask));
%     stdQ(ii) = std(Q(mask));
%     
% end
%%
% figure;
% plot(t_vec2(1:ii-1),2*meanQ(1:ii-1)./stdQ(1:ii-1).^2)
% xlabel('t');
% ylabel('2<Q>/ var(Q)');
% figure;
% plot(t_vec2(1:ii-1),meanQ(1:ii-1))
% xlabel('t');
% ylabel('<Q>');
% 
% figure;
% plot(t_vec2(1:ii-1),stdQ(1:ii-1).^2)
% xlabel('t');
% ylabel('Var(Q)');
% 
%%
%calculate mean and sig of Q within time steps of size dt
% dt=diff(t_vec,1,1);
% dt = mean(dt(:));
% t=0;
% ii=1;
% while t<t_min
%     mask = [abs(diff(t_vec<t,1,1));zeros(1,numRuns)]~=0;
%
%     meanQ(ii) = mean(Q(mask));
%     stdQ(ii) = std(Q(mask));
%     t = t+dt;
%   ii=ii+1;
% end
% %%
% figure;
% plot((1:ii-1)*dt,meanQ./((1:ii-1)*dt));
% figure;
% plot((1:ii-1)*dt,stdQ./((1:ii-1)*dt));
%
%
% figure;
% plot((1:ii-1)*dt,2*meanQ(1:ii-1)./stdQ(1:ii-1).^2);
%%
