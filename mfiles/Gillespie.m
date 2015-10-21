% function P = Gillespie(w_p,w_m)
clear site_t t_vec
N=length(w_p);
numRuns = 1e3;

Q=zeros(1,numRuns);
%
% W = -diag(w_p + w_m([N,1:N-1]),0) + ...
%     diag(w_m(1:N-1),1) + ...
%     diag(w_p(1:N-1),-1);
% W(N,1) = w_m(N);
% W(1,N) = w_p(N);
% lambda = sort(real(eig(W)));
% lambda_max = lambda(end-1);
% t_end = -1/lambda_max;
t_end=1/min([w_p,w_m]);
r = w_p + w_m([N,1:N-1]);
for iR=1:numRuns
    wb=waitbar(iR/numRuns);
    
    t=0;
    site = 0;
    iT=1;
%     tic;
    while (t<t_end*10000)
        n = mod(site,N)+1;
        n_1 = mod(site-1,N)+1;
        
        a = rand(2,1);
        
        delta_t = -1/r(n)*log(a(1)/r(n));
%         delta_t_vec(iT,iR) = delta_t;
        
        opts = [w_p(n),w_m(n_1)];
        [w1,I1] = min(opts);
        [w2, I2] =  max(opts);
       
        if(a(2)<w_p(n)/r(n))
            site=site+1;
        else
            site=site-1;
        end
        
%         if (a(2)<w1/r(n))
%            
%             site = site - (2*I1-3);
%             
%         else
%             
%             site = site - (2*I2-3);
%           
%         end
%         
        site_t(iT,iR) = site;%mod(site,N)+1;
        
        t = t + delta_t;
        t_vec(iT) = t;
        iT = iT+1;
        
        
        
        
        
    end
%         toc;
% 
% 
%     ind1 = find((mod(site_t(:,iR),N)+1)==4);
%     if(~isempty(ind1) && ind1(end)==iT-1)
%         ind1 = ind1(1:end-1);
%     end
%     if(~isempty(ind1) && ind1(1)==1)
%         ind1 = ind1(2:end);
%     end
%     ind2 = find(mod(site_t(ind1-1,iR),N)+1==3);
%     ind3 = find(mod(site_t(ind1+1,iR),N)+1==3);
%     
%     Q(iR) = length(ind2)-length(ind3);
%     muQ(iR,iT)

end
%%
idx=1;
for iT = 2:length(t_vec)
    [mu(idx),sig(idx),Q(idx,:)] = Q_stats (site_t,1,iT,N);
    idx=idx+1;

end
%%
figure;plot(t_vec(2:iT-1),2*mu(1:idx-1)./sig(1:idx-1).^2);
x=-max(abs(Q(end,:))):1:max(abs(Q(end,:)));

P1 = hist(Q(end,:),x);
P2 = hist(-Q(end,:),x);
figure;plot(x,P1,x,normpdf(x,mu(end),sig(end))*numRuns)
figure;plot(x,log(P1./P2))

%%

% %%
% figure;
% plot(sort(Q),1:length(Q))
% mu = mean(Q);
% sig= std(Q);
% x = linspace(min(Q),max(Q),100);
% hold on
% plot(x, 0.5*(1+erf((x-mu)/sig/sqrt(2)))*length(Q))
% 
% figure;
% plot(sort(Q),length(Q):-1:1)
% mu = mean(Q);
% sig= std(Q);
% x = linspace(-max(abs(Q)),max(abs(Q)),100);
% hold on
% plot(x, 0.5*(1-erf((x-mu)/sig/sqrt(2)))*length(Q))
% 
% 
figure;
mu1=mu(end);
sig1=sig(end);
[h,x]=hist(Q(end,:),-max(abs(Q(end,:))):1:max(abs(Q(end,:))));
hist(Q(end,:),-max(abs(Q(end,:))):1:max(abs(Q(end,:))));
hold on
plot(x,1/sqrt(2*pi*sig1^2)*exp(-(x-mu1).^2/(2*sig1^2))*length(Q(end,:)),'r','LineWidth',4)
plot(-max(abs(Q(end,:))):1:max(abs(Q(end,:))),h)
% end