% clear all
sigma =5;
Delta=1;
N=10;
beta=0.1;
w_beta=1;

SMF_vec=[linspace(0,5,10), linspace(5,85,10)];
% SMF_vec=20;
% SMF_vec=0:5;
for q=1:length(SMF_vec)
    
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
%     w_p = exp(2*SMF/((N-1))*(0:N-1));
    
%     
%     w_p = [3,20,3,20,3,20,3,20,3,20]
% %     w_p = [exp(5*SMF/4),exp(-1/4*SMF),exp(5*SMF/4),exp(-1/4*SMF),exp(5*SMF/4),exp(-1/4*SMF),exp(5*SMF/4),exp(-1/4*SMF),exp(5*SMF/4),exp(-1/4*SMF)]
%     w_p = [exp(SMF/3),exp(2/3*SMF),exp(SMF/3),exp(2/3*SMF),exp(SMF/3),exp(2/3*SMF),exp(SMF/3),exp(2/3*SMF),exp(SMF/3),exp(2/3*SMF)]
%     w_m = ones(1,N);
%     Gillespie3(w_p,w_m);
%     
% %     w_m= [2,4,2,4,2,4,2,4,2,4]
%     [muQ(q),sigQ(q),muX(q),sigX(q),t_end(q)] = Gillespie3(w_p,w_m);
% figure;plot(cumsum(log(w_p./w_m)))%     Delta1 = 1;  
% as(q) = (w_p(1)+w_m(2))/(w_m(1)+w_p(2));
% mu(q) = 2*(w_p(1)*w_p(2)-w_m(1)*w_m(2))/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
% D(q) = (2*(w_p(1)*w_p(2)+w_m(1)*w_m(2)) - mu(q)^2)/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
% 

 w_p = [exp(3/2*SMF),exp(-SMF/2),exp(3/2*SMF),exp(-SMF/2),exp(3/2*SMF),exp(-SMF/2),exp(3/2*SMF),exp(-SMF/2),exp(3/2*SMF),exp(-SMF/2)];
    w_m = ones(1,N);
    [muQ3(q),sigQ3(q),muX3(q),sigX3(q),t_end3(q)] = Gillespie3(w_p,w_m);
% figure;plot(cumsum(log(w_p./w_m)))%     Delta1 = 1;  
as3(q) = (w_p(1)+w_m(2))/(w_m(1)+w_p(2));
mu3(q) = 2*(w_p(1)*w_p(2)-w_m(1)*w_m(2))/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
D3(q) = (2*(w_p(1)*w_p(2)+w_m(1)*w_m(2)) - mu3(q)^2)/(w_p(1)+w_p(2)+w_m(1)+w_m(2));


%     A1=exp(SMF)
%     B1=1
%     
%     w_p=[A1,B1,A1,B1,A1,B1,A1,B1,A1,B1];
%     w_m=1./w_p;
%     [muQ(q),sigQ(q),t_end(q)] = Gillespie3(w_p,w_m);
% %   
%     r = w_p + w_m([N,1:N-1]);
%     t_end(q) = max(1./r)*1e2;


%     w_p=ones(1,N);
%     w_m=w_p;
%     [muQ(q),sigQ(q),t_end(q)] = Gillespie3(w_p,w_m);
%     
%     
%     w_p=[1,100,1,100,1,100,1,100,1,100]
%     w_m=w_p;
%     [muQ2(q),sigQ2(q),t_end2(q)] = Gillespie3(w_p,w_m);
%     
    
    %bimodal chain
    w_p = ones(1,N)*exp(SMF/2);
    w_m = ones(1,N)*exp(-SMF/2);
% figure;plot(cumsum(log(w_p./w_m)))
SMF2(q) = sum(log(w_p./w_m))/N;
    [muQ2(q),sigQ2(q),muX2(q),sigX2(q),t_end2(q)] = Gillespie3(w_p,w_m);
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
SMF=0:85;

plot(SMF, 2*tanh(SMF/2),'--b','LineWidth',3);

plot(SMF, 2*N*tanh(SMF/2/N)./(1+(tanh(SMF/2/N).^2+tanh(SMF/2/N).^3+...
+tanh(SMF/2/N).^4+tanh(SMF/2/N).^5+tanh(SMF/2/N).^6+...
    +tanh(SMF/2/N).^7+tanh(SMF/2/N).^8+tanh(SMF/2/N).^9+tanh(SMF/2/N).^10)),'r','LineWidth',3)


plot(SMF2,2*muQ2./sigQ2.^2/N,'b*','MarkerSize',10); %1 site per cell
% plot(abs(SMF2),2*abs(muQ3)./sigQ3.^2,'ro','LineWidth',4); %10 sites per cell


xlabel('s');
ylabel('\mu_M / D_M');
legend(['M=1 ';...
        'M=10'],   'Location', 'SouthEast')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/mu_sig_ratio2.eps');
% print(gcf, '-depsc2', '/users/physics/hurowits/MyWebSite/PROJ/NEK/Figs/mu_sig_ratio.eps');

%%
figure;
axes('FontSize',24);
hold on;
grid on
% plot(SMF2,muQ2./sigQ2.^2,'*')

plot(SMF2,muX2./sigX2.^2,'b*','MarkerSize',10); %bi modal chain
plot(abs(SMF_vec),abs(muX)./sigX.^2,'ro','LineWidth',3); %random chain

plot(0:30,tanh((0:30)/2),'--b','LineWidth',3);

plot(abs(SMF_vec),abs(mu)./(D*2),'r','LineWidth',3); %random chain


plot(SMF_vec, 2*tanh(SMF_vec/4)./(1+tanh(SMF_vec/4).^2))

% plot(abs(SMF_vec),abs(mu3)./(D3*2),'--g','LineWidth',3); %random chain
% plot(abs(SMF_vec),SMF_vec/2,'k','LineWidth',3); %random chain


xlabel('s');
ylabel('<x> / var(x)');
legend(['asymetric & invariant ';...
        'asymetric & dispersive';,...
        
        'tanh (s/2) [R-N]      ';...
        '\mu/2D                '],...
    'Location', 'SouthEast');
axis tight;
% axis([0 30,0, 1.2])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/mu_sig_ratio3.eps');