clear all
sigma =5;
Delta=1;
N=10;
beta=0.1;
w_beta=1;

SMF_vec=[0:5, linspace(6,100,5)];
% SMF_vec=20;
for q=1:length(SMF_vec)
%     q=7
    SMF = SMF_vec(q);
    
    %     while abs(SMF)<=SMF_vec(q)
    
    %         Energy = randn(1,N)*Delta;
    %         Energy=1:N;
    %         G0 = linspace(0,1,N);
    %         randind = randperm(N);
    %         G0 = G0(randind);
    %         G = 10.^(G0*sigma - sigma/2);
    %         G = G*mean(1./G);
    %
    % engineer transition rates for increased SMF
    % for p=1:N-1
    %     if (log10(G(p)*0.01)<0)
    %
    %         E1 = max(Energy(p:p+1));
    %         E2 = min(Energy(p:p+1));
    %         Energy(p)=E1;
    %         Energy(p+1)=E2;
    %     end
    % end
    %
    %
    %         w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
    %         w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];
    %
    %         w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
    %         w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];
    %
    
    
    %
    %         w_drv =1e-5*G;
    % w_drv =0.2719*G;
    %         % w_drv =300*G;
    %         w_drv=0.01*G;
    %         w_m = (w_drv + w_bath_ccw); %w_{n-1,n}
    %         w_p = (w_drv + w_bath_cw);  %w_{n,n-1}
    
    %"random" chain
%     w_p = exp(SMF/(N*(N-1))*(0:N-1));
%     w_m = exp(-SMF/(N*(N-1))*(0:N-1));
%     [muQ(q),sigQ(q),t_end(q)] = Gillespie3(w_p,w_m);
%     
%     
   
%     r = w_p + w_m([N,1:N-1]);
%     t_end(q) = max(1./r)*1e2;


    w_p=[5,30,5,30,5,30,5,30,5,30];
    w_m=[1,2,1,2,1,2,1,2,1,2];
    
    A1=exp(Delta1)
    B1=exp(SMF-Delta1)
    
    w_p=[A1,B1,A1,B1,A1,B1,A1,B1,A1,B1];
    w_m=1./w_p;
%     
%     w_p=[exp(g),exp(d),exp(g),exp(d),exp(g),exp(d),exp(g),exp(d),exp(g),exp(d)]
%     w_m=1./w_p;
%     
    mu=(w_p(1)*w_p(2)-w_m(1)*w_m(2))/(w_m(1)+w_m(2)+w_p(1)+w_p(2))
    sig =sqrt(((w_p(1)*w_p(2)+w_m(1)*w_m(2))-2*mu^2)/(w_m(1)+w_m(2)+w_p(1)+w_p(2)))
    
    [muQ(q),sigQ(q),t_end(q)] = Gillespie3(w_p,w_m);
    
    w_p=ones(1,N);
    w_m=w_p;
    [muQ(q),sigQ(q),t_end(q)] = Gillespie3(w_p,w_m);
    
    
    w_p=[1,100,1,100,1,100,1,100,1,100]
    w_m=w_p;
    [muQ2(q),sigQ2(q),t_end2(q)] = Gillespie3(w_p,w_m);
    
    w_p = 1:N;
    w_m = w_p;
        [muQ3(q),sigQ3(q),t_end3(q)] = Gillespie3(w_p,w_m);

    %bimodal chain
%     w_p = ones(1,N)*exp(SMF/2);
%     w_m = ones(1,N)*exp(-SMF/2);
%     SMF2(q) = sum(log(w_p./w_m))/N;
%     [muQ2(q),sigQ2(q)] = Gillespie3(w_p,w_m);
%     
    SMF
    q
end
%%

figure;
axes('FontSize',24);
hold on;
grid on
% plot(SMF2,muQ2./sigQ2.^2,'*')
plot(abs(SMF_vec),abs(muQ)./sigQ.^2,'ro','LineWidth',2); %random chain

plot(SMF2,muQ2./sigQ2.^2/N,'b*'); %bi modal chain

plot(0:100,tanh((0:100)/2),'g','LineWidth',2);
xlabel('SMF');
ylabel('\mu_Z / \sigma_Z^2N');
legend(['Random chain  ';...
    'Bi-modal chain';...
    'tanh (SMF/2)  '],...
    'Location', 'SouthEast')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/mu_sig_ratio.eps');
% print(gcf, '-depsc2', '/users/physics/hurowits/MyWebSite/PROJ/NEK/Figs/mu_sig_ratio.eps');

