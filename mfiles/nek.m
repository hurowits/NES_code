% Script for calculating current and statistics for nearest neighbor jump
% process on a ring.

clear all

% -----------------------
% Running parameters
% -----------------------

numRealizations=1;
flagAnalytic=0;  %0 for numerical solution, 1 for analytic solution

% -----------------------
% Model parameters
% -----------------------

N=100; %number of sites

Delta = 1;                             %Energy interval
w_beta = 1;                            %Coupling to bath
beta = 0.1;                            %Inverse temperature
sigma_vec =3;                        %log-width of distribution


nu_vec= logspace(-8,5,100);  %Driving intensity
harmonicMeanG = zeros(length(sigma_vec),1);

%-----------------------
% Variable initialization
%-----------------------

w_N0 = zeros(1,N);
w_0N = zeros(1,N);
P = zeros(1,N);
A = zeros(1,N);
P_inf=1/N*ones(1,N);
for iR = 1:numRealizations
    wb=waitbar(iR/numRealizations);
    
    Energy = randn(1,N)*Delta;
    P0 = exp(-Energy*beta)/sum( exp(-Energy*beta));
%     G0 = rand(1,N);
    G0 = linspace(0,1,N);
    randind = randperm(N);
    G0 = G0(randind);
    for iS = 1:length(sigma_vec)
        
        sigma = sigma_vec(iS);
        G = 10.^(G0*sigma - sigma/2);
%         G = rand(1,N);
        G = G*mean(1./G);
        
        Delta_n = Energy(2:N)-Energy(1:N-1);
        Delta_n = [Delta_n, Energy(1)-Energy(N)];
        
        
        w_bath_cw = 2*w_beta./(1+exp(-beta*(Energy(1:N-1)-Energy(2:N)))); %w_plus
        w_bath_cw = [w_bath_cw  2*w_beta/(1+exp(-beta*(Energy(N)-Energy(1))))];
        
        w_bath_ccw = 2*w_beta./(1+exp(beta*(Energy(1:N-1)-Energy(2:N)))); %w_minus
        w_bath_ccw = [w_bath_ccw,  2*w_beta/(1+exp(beta*(Energy(N)-Energy(1))))];
        
        
        for iE = 1: length(nu_vec)
            
            w_drv = nu_vec(iE) * G * w_beta;
            w_m = (w_drv + w_bath_ccw); %w_{n-1,n}
            w_p = (w_drv + w_bath_cw);  %w_{n,n-1}
%             
            
               
                W_cl = -diag(w_drv + w_bath_cw +w_drv([N,1:N-1])+ w_bath_ccw([N,1:N-1]),0)+...
                    diag(w_drv(1:N-1) + w_bath_ccw(1:N-1),1) + ...
                    diag(w_drv(1:N-1) + w_bath_cw(1:N-1),-1);
                
                W_cl(N,1) = w_bath_ccw(N) + w_drv(N);
                W_cl(1,N) = w_bath_cw(N) + w_drv(N);
                
                I_mat = -diag(w_m([N,1:N-1]) - w_p,0)/N;
                
                
                B = [W_cl;ones(1,N)];
                C = zeros(N+1,1);C(N+1)=1;
                [P_cl,R] = linsolve(B,C);
%               
                
            [P_an1,P_an2,I1(iE),er1(iE)] = ness(w_m,w_p);
            er(iE) = sum(abs(P_an1-P_an2));
            [P1,I1(iE),er1(iE)] = ness(w_m,w_p,0);
            
%             er3(iE)=max(W_cl*P1');
%             er3(iE) = sum((P1-P_cl').^2);

            V = [0 cumsum(log(w_m./w_p))];
            Vflat = V-V(end)*(0:N)/(N);
            Vflat1 =  Vflat(1:end-1);
            
            I_exact(iE) = P1(1)*w_p(1)-P1(2)*w_m(1);         
            
            SMF(iE)=-V(end);
            
            D_R_1(iE) =P_inf* I_mat *pinv(W_cl)*I_mat *P1';
            I_1(iE) = P_inf*I_mat*P1';
            
            I_fdt1(iE) = 0.5*(w_p(1)*P1(1)+w_m(1)*P1(2) - I1(iE)^2)*(SMF(iE))/beta;
            
            D_p(iE) = 0.5*((w_p*P1'+w_m*P1([2:N,1])')/N^2 - I_exact(iE)^2);
            D_p0(iE) = 0.5*((w_bath_cw*P1'+w_bath_ccw*P1([2:N,1])')/N^2 - I_exact(iE)^2);
%             D_p0_a(iE) = 0.5*((w_bath_cw(1)*P1(1)+w_bath_ccw(1)*P1(2)) - I_exact(iE)^2);

            %             if iE==1
            %                 P0=P1;
            %             end
            % do again for nu+dnu
            dnu = 1e-5;
            w_drv = (nu_vec(iE)+dnu) * G * w_beta;
            w_m = (w_drv + w_bath_ccw); %w_{n-1,n}
            w_p = (w_drv + w_bath_cw);  %w_{n,n-1}
            [P2,I2(iE),er2(iE)] = ness(w_m,w_p,0);
            
            V = [0 cumsum(log(w_m./w_p))];
            Vflat = V-V(end)*(0:N)/(N);
            Vflat2 =  Vflat(1:end-1);
            SMF2(iE) = -V(end);
            
            BG = -(Vflat2-Vflat1)/dnu;
            
            I_fdt2(iE) = 0.5*(w_p(1)*P2(1)+w_m(1)*P2(2) - I2(iE)^2)*(SMF2(iE))/beta;
            chi_fdt(iE) = (I_fdt2(iE)-I_fdt1(iE))/dnu;
            
            B = (log(P2)-log(P1))/dnu;
            B2 = 1./P1.*(P2-P1)/dnu;
%             B_lrt1 = (sum(Delta.*G)-N*G(1)*Delta_n(1))*beta/sum(exp(beta*Energy))
             
            temp(iE) = B(1)-B(2);
            
            chi_numerical(iE) = (I2(iE)-I1(iE))/dnu;
            
            chi(iE) = G(1)*(P1(1)-P1(2))+...
                w_p(1)*(P2(1)-P1(1))/dnu - w_m(1)*(P2(2)-P1(2))/dnu ;
            
            chiG(iE) = G(1)*(P1(1)-P1(2))+...
                w_p(1)*BG(1)*P1(1) - w_m(1)*BG(2)*P1(2);
            
            Q_dot1(iE) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P1'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P1([N,1:N-1])' ; %equation (31)
            
            Q_dot2(iE) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*P2'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*P2([N,1:N-1])' ; %equation (31)
            
            chi_Q_numerical(iE) = (Q_dot2(iE)-Q_dot1(iE))/dnu;
            
            chi_Q(iE) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*(P1.*B)'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*(P1([N,1:N-1]).*B([N,1:N-1]))';
            
            
            chi_Q_G(iE) = (w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1]))*(P1.*BG)'...
                - (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1]))*(P1([N,1:N-1]).*BG([N,1:N-1]))';
            
        end
    end
end

%%
% 
% D_R = I1./SMF/2;
% 
% figure;
% axes('FontSize',16);
% hold on
% grid on
% plot(nu_vec,D_R,nu_vec,abs(D_R_1),nu_vec,abs(D_R_1)./D_R,nu_vec,abs(I_1),nu_vec,abs(SMF),...
%     nu_vec,D_p0,'LineWidth',4);
% plot([1/max(G) 1/max(G)], [1e-8 1e8],'--k')
% plot([1/min(G) 1/min(G)], [1e-8 1e8],'--k')
% 
% xlabel('nu');
% legend(['D_R FDT ';...
%         'D_R corr';...
%         'ratio   ';...
%         'abs(I)  ';...
%         'abs(SMF)']);
% 
% axtype(3)
%%
chi_eq = P0(2)*w_bath_cw(2)*((sum(G.*Delta_n)-N*Delta_n(2)*G(2))/sum( exp(-Energy*beta)) +Delta_n(2)*G(2))*beta
chi_eq = P0(1)*w_bath_cw(1)*((sum(G.*Delta_n)-N*Delta_n(1)*G(1))/sum( exp(-Energy*beta)) +Delta_n(1)*G(1))*beta
chi_eq=0;
for ii=1:N
chi_eq =chi_eq+ P0(ii)*w_bath_cw(ii)*((sum(G.*Delta_n)-N*Delta_n(ii)*G(ii))/sum( exp(-Energy*beta)) +Delta_n(ii)*G(ii))*beta;
end
chi_eq=chi_eq/N
% chi_eq = 1/N*((sum(G.*Delta_n)-N*Delta_n(1)*G(1))/N+Delta_n(1)*G(1))*beta

chi_0 = sum(G.*Delta_n)*beta/N^2;
chi_inf =  sum(Delta_n./G.*P2)*beta/N^2;


chi_0_2 =SMF./nu_vec.*D_p;
% chi_0_3 =SMF./nu_vec.*D_p0_a;

chi_Q_2 = 0.5* ((w_bath_ccw([N,1:N-1]).*Delta_n([N,1:N-1])).^2*P0'...
                + (w_bath_cw([N,1:N-1]).*Delta_n([N,1:N-1])).^2*P0([N,1:N-1])' ); 
% chi_0_2 = 0.5*beta/N* (w_bath_ccw([N,1:N-1])*P0'...
%                 + w_bath_cw([N,1:N-1])*P0([N,1:N-1])' )*V(end)/nu_vec(1); 
%              
% chi_0 = (G(1)*(P0(1)-P0(2))+B(1)*w_bath_ccw(1)*P0(1)-B(2)*w_bath_cw(1)*P0(2))/N;


figure;
axes('FontSize',22);
hold on;
grid on
plot(nu_vec,chi,'g*','MarkerSize',10)

% plot(nu_vec,X,nu_vec,Y,nu_vec,X+Y)
plot(nu_vec,chi_numerical,'ro','MarkerSize',10)

plot(nu_vec,chiG,'bs','MarkerSize',10)

% plot(nu_vec,I_kubo,'g');
% plot(nu_vec,chi_0_2,'dc','MarkerSize',10)
plot(nu_vec,chi_fdt/N*beta,'dc','MarkerSize',10)
% plot(nu_vec,I_exact .* abs((chi_inf-chi_0)./abs(max(I_exact)-min(I_exact))),'--k','LineWidth',4);
% plot(nu_vec,I_fdt1 .* abs((chi_inf-chi_0)./abs(max(I_fdt1)-min(I_fdt1))),'--c','LineWidth',4);

plot(nu_vec,chi_eq*ones(1,length(nu_vec)),':r','LineWidth',4)
plot(nu_vec,chi_0*ones(1,length(nu_vec)),'--r','LineWidth',4)

% plot(nu_vec,SMF/N^2./nu_vec,':g','LineWidth',4)

% plot(nu_vec,chi_inf*ones(1,length(nu_vec)),'--b','LineWidth',4)
% plot(nu_vec, beta/N^2*sum(1./G.*Delta_n)*ones(1,length(nu_vec)),'--r','LineWidth',4)

% plot(nu_vec(2:end),diff(I1)./diff(nu_vec))
xlabel('\nu');
ylabel('\chi_Z  ')
legend(['Analytical        ' ;...
        'Numerical         ';...
        'canonical-like FDT';
        'canonical FDT     '],'Location','SouthWest');
       % 'Current           ';...
        %'Current_{fdt}     ']
        
%     axis([nu_vec(1) nu_vec(end) 2*min([chi(:)',chi_fdt(:)',chiG(:)']) 2*max([chi(:)',chi_fdt(:)',chiG(:)'])])
v=axis;
axis([nu_vec(1), nu_vec(end) v(3) v(4)])
axtype(1)

% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/chi_Z.eps');

%%

chi_Q_0 = sum(P0.*G.*Delta_n.^2)*beta;
chi_Q_inf = (((mean(Delta_n.^2./G))*beta./nu_vec.^2));
figure;
axes('FontSize',24);
hold on;
grid on
plot(nu_vec,chi_Q,'g*')

% plot(nu_vec,X,nu_vec,Y,nu_vec,X+Y)
plot(nu_vec,chi_Q_numerical,'ro')
% plot(nu_vec,I_kubo,'g');

% plot(nu_vec,chi_Q_2(1)*ones(1,length(nu_vec)),'-r','LineWidth',4)
plot(nu_vec,chi_Q_G,'sb','MarkerSize',10)
plot(nu_vec,Q_dot1,'--k','LineWidth',4);
plot(nu_vec,chi_Q_0*ones(1,length(nu_vec)),'--r','LineWidth',4)
plot(nu_vec,chi_Q_inf,'--b','LineWidth',4)

% plot(nu_vec,chi_inf*ones(1,length(nu_vec)),'--b','LineWidth',4)
% plot(nu_vec, beta/N^2*sum(1./G.*Delta_n)*ones(1,length(nu_vec)),'--r','LineWidth',4)

% plot(nu_vec(2:end),diff(I1)./diff(nu_vec))
xlabel('\nu');
ylabel('\chi_Q, dQ/dt')
legend(['Analytical  ' ;...
        'Numerical   ';...
        'Non-eq FDT  ';
        'Heat current'],'Location','SouthWest');
axtype(3)
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEK/Figs/chi_Q.eps');