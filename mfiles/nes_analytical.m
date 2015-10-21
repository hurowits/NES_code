clear all

cols = 'rgbckym';c=0;
M=[10];

a_0=1;L=a_0*M;
for iR=1:1

        c=c+1;
        a=1;b=1;
        G = ones(1,M);
        % s0 = rand(1,M)*2-1;
        %  s0 = s0/sum(s0);3
        sigma=1;
        s0 = (linspace(0,1,M)*2-1);
        %     s0 = rand(1,M)*2-1;
        
        SMF_vec = [linspace(0,0.002,100),linspace(0.2,30*M,1000)];
%         SMF_vec=linspace(-30*M,30*M,1000);
%         SMF_vec=1000;
        %     SMF_vec=500;
        % s0 =2*(-M/2:M/2)/(M/2*(M/2-1));
        randInd1 = randperm(M);
        randInd2 = randperm(M);
        % s0=s0(randInd);
        
        % SMF_vec=0;
        % s0=1;
        sigma_vec=[linspace(0,10,5)];
        sigma_vec=[0,2.8];
%         sigma_vec=[0,6];
%         sigma_vec=linspace(0,20,15);
        %     sigma_vec=5;
        for iSig = 1:length(sigma_vec);
            sigma = sigma_vec(iSig);
            for q=1:length(SMF_vec)
                SMF=SMF_vec(q);
                
                %         s = sigma*s0+SMF/(M*(M-1))*(0:M-1)*2;
                %         s_p = sigma*s0(randInd1)+SMF/(M*(M-1))*(0:M-1)*2;
                %         s_m = -(sigma*s0(randInd2)+SMF/(M*(M-1))*(0:M-1)*2);
                s_p = sigma*s0(randInd1)+SMF/M*ones(1,M);
                s_m = -(sigma*s0(randInd2)+SMF/M*ones(1,M));
                
                
                %         %     w_p = exp(-SMF/(M*(M-1))*(0:M-1));
                %     w_m = exp(SMF/(M*(M-1))*(0:M-1));
                %     %     s1 =
                
                w_p = G.*exp(s_p/2);
                w_m = G.*exp(-s_p/2);
                
                %                     figure(1);
                %                     hold on
                %                     plot(cumsum(log(w_p./w_m)));
                %                     title(num2str(SMF_vec(q)));
                %                     axis([1 M -700 700])
                
                r = w_p + w_m([M,1:M-1]);
                t_end = max(1./r);
                t=t_end;
                
                
                W = -diag(w_p + w_m([M,1:M-1]),0) + ...
                    diag(w_m(1:M-1),1)+diag(w_p(1:M-1),-1);
                W(1,M) = w_p(M);
                W(M,1) = w_m(M);
                
                
                
                delta_k = 0.1;
                
                if (M==2)
                    
                    W(1,2) = w_p(2) + w_m(1);
                    W(2,1) = w_m(2) + w_p(1);
                    [g0,lambda0,C0] = charFn(W,0);
                    
                    W(1,2) = w_p(2)*exp(1i*delta_k)  + w_m(1);
                    W(2,1) = w_m(2)*exp(-1i*delta_k) + w_p(1);
                    [g1,lambda1,C1] = charFn(W,0);
                    
                    W(1,2) = w_p(2)*exp(1i*2*delta_k)  + w_m(1);
                    W(2,1) = w_m(2)*exp(-1i*2*delta_k) + w_p(1);
                    [g2,lambda2,C2] = charFn(W,0);
                else
                    
                    [g0,lambda0,C0] = charFn(W,0);
                    [g1,lambda1,C1] = charFn(W,delta_k);
                    [g2,lambda2,C2] = charFn(W,2*delta_k);
                    
                end
                
                %     dk_vec=linspace(0,3,10000);
                %     for k=1:length(dk_vec)
                %         dk=dk_vec(k);
                %         [g_dbg(k),lambda_dbg(k),C_dbg(k)] = charFn(W,dk);
                %
                %         if (M==2)
                %
                %             W(1,2) = w_p(2)*exp(1i*dk)+w_m(1);
                %             W(2,1) = w_m(2)*exp(-1i*dk) + w_p(1);
                %             [g_dbg(k),lambda_dbg(k),C_dbg(k)] = charFn(W,0);
                %
                %         end
                %     end
                %     %     figure(1);
                %     hold on;
                %     plot(dk_vec,lambda_dbg,dk_vec,C_dbg,dk_vec(1:end-1),diff(C_dbg).^2./dk_vec(1:end-1).^2,...
                %         dk_vec(1:end-2),diff(diff(C_dbg))./dk_vec(1:end-2).^2,...
                %         dk_vec(1:end-1),diff(lambda_dbg)./dk_vec(1:end-1));
                %
                
                %
                dCdk(q) = (C1-C0)/delta_k;
                dlambdadk = (lambda1 - lambda0)/delta_k;
                %
                % %     dlambdadk = (lambda_dbg(2) - lambda_dbg(1))/(dk_vec(2)-dk_vec(1));
                %     mobility(q) = -1i*dlambdadk;
                %     Diffusion(q) = -0.5*(lambda2 - 2*lambda1 + lambda0)/delta_k^2 - dCdk(q)*dlambdadk;
                %
                mobility(q,iSig) = -i*dlambdadk;
                Diffusion(q,iSig) = -0.5*(lambda2 - 2*lambda1 + lambda0)/delta_k^2 ;
                
                v_drda(q,iSig) = (1- mean(w_m./w_p))./(mean(1./w_p));
                D_drda(q,iSig) = 0.5*(mean(1./w_p)).^(-3)*...
                    (1-(mean(w_m./w_p)).^2)./...
                    (1-mean((w_m./w_p).^2))*...
                    (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
                    2*mean(w_m./w_p.^2)*mean(1./w_p));
                
                
                
                %
                %     sigma(q) = 1- 4*a*b * cosh((s1+s2)/4).^2./(a*cosh(s1/2)+b*cosh(s2/2)).^2;
                %     rat2(q) = 4*tanh(SMF/4)/(1+sigma(q) * tanh(SMF/4)^2);
                %
                %
                
                %     varX(q) = -(g2 - 2*g1 + g0)/delta_k^2- meanX(q)^2;
                
                %     meanX(q) = -1i*(g1 - g0)/delta_k;
                %     varX(q) = -(g2 - 2*g1 + g0)/delta_k^2- meanX(q)^2;
            end %for SMF
            a_inf = mean((1./w_p).^2)/mean(1./w_p)^2;
        end %for sigma
        f(iR,:)=-(real(mobility(:,2))./real(Diffusion(:,2)))/2;
        
end
%%
% figure;plot(SMF_vec,(meanX));
% figure;plot(SMF_vec,abs(varX));
% figure;plot(SMF_vec,2*(meanX)./(varX))
s = SMF_vec/M;
sigma=sqrt(sigma_vec(2));
V = 2*exp(sigma/2)*sinh(s-sigma);
D=exp(sigma)*sinh(s-sigma)./sinh(s-2*sigma).*...
(0.5*exp(s+sigma/2) + 0.5*exp(-s).*(2*exp(7/2*sigma)-exp(5/2*sigma)));

figure;
axes('FontSize',24);
grid on
hold on;
plot(s, v_drda./D_drda,'--','LineWidth',4);
plot(s, -real(mobility./Diffusion)/M,'LineWidth',2)

%%
SMF_mat = repmat(SMF_vec,iR,1);

a_s3 = 1+(a_inf-1)*tanh(s).^2;
a_s =( 2./(1+exp(s))+a_inf./(1+exp(-s)));
a_s2 =( 2./(1+exp(s/2))+a_inf./(1+exp(-s/2)));
a_s4 = (1+s*a_inf)./(1+s);
figure;
% 
axes('FontSize',24);
grid on
hold on;
plot(s,s/2,'b--','LineWidth',4)

plot(s,-real(mobility(:,1))./real(Diffusion(:,1))/M/2,'k','LineWidth',4);
plot(s,1/M*tanh(SMF_vec/2),'k','LineWidth',4);

plot(s,f(5:15,:)/M,'r-','LineWidth',1);
% plot(s,1./(2*f/M./tanh(SMF_mat/2/M)),'g--','LineWidth',1)
% plot(s,2*f/M)
% plot((SMF_vec/M),f([5:15],:)/M,'r-','LineWidth',1)
% plot((SMF_vec/M),mean(f)/M,':r','LineWidth',4)
plot(s,1/a_inf * tanh(a_inf* s/2),'k','LineWidth',4)
% plot(s, 2./a_s.*tanh(s/2),'k','LineWidth',4);
% plot(s, 2./a_s2.*tanh(s/2),'--k','LineWidth',4);
% plot(s, 2./a_s3.*tanh(s/2),'--r','LineWidth',4);
% plot(s, 1./(2./a_s4.*tanh(s/2)),'--g','LineWidth',4);
% a_inf_2=sigma*coth(sigma);
% plot((SMF_vec/M),2/a_inf_2 * tanh(SMF_vec/2),'c','LineWidth',4)

% legend(['\sigma=0     ';...
%         '\sigma=\infty']);
% legend(num2str(sigma_vec'));
% legend(['ESR';...
%         '\sigma = 0     ';...
%         '\sigma = \infty';...
%         '\
xlabel('s');
ylabel('v/2D');
axis([10^-2 2*M, 0 1.1]);
axtype(1);

% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/v_D_rat.eps');
%%
figure;
plot(s, tanh(s/2)./(2*f(1,:)/M))
%%
%     figure;
%     axes('FontSize',24);
%     grid on
%     hold on;
% %     plot(tanh(SMF_vec/2/M_vec(1)),-(real(Diffusion)./real(mobility)).*repmat(tanh(SMF_vec/2/M),length(sigma_vec),1)');
%     plot(tanh(SMF_vec/2/M_vec(1)),squeeze(-(real(Diffusion(1,:,1:15))./real(mobility(1,:,1:15)))).*repmat(tanh(SMF_vec/2/M_vec(1)),15,1)','--');
%     plot(tanh(SMF_vec/2/M_vec(1)),tanh(SMF_vec/2/M_vec(1)),'--k','LineWidth',2)
%     plot(tanh(SMF_vec/2/M_vec(1)),1/M_vec(1)*ones(1,length(SMF_vec)),'-k','LineWidth',2)
%
% %     legend([num2str(sigma_vec');...
% %         '2Mx     ']);
%     xlabel('tanh(s/2M)');
%     ylabel('tanh(s/2M) D/\mu ');
%
%     xlabel('x');
%     ylabel('x D/\mu ');
%
%     axis([0 1.0, 0 1.0])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/x_D_mu_rat.eps');

%%
figure;
% hold on
axes('FontSize',24);
grid on;
hold on;
plot(sigma_vec,M./(-real(mobility(end,:))./real(Diffusion(end,:))/2),'ro','MarkerSize',10,'LineWidth',4);

s=sigma_vec;
% s=linspace(0,20,1000);

for iS=1:length(s)
    sigma=sigma_vec(iS);
    w=exp(sigma*s0(randInd1)+SMF/M*ones(1,M));
    h1=1/mean(1./sqrt(w));
    h2=1/(mean(1./w));
    F(iS) = 1/( h1^2/h2);
end


plot(sigma_vec,1./F,'.-b','LineWidth',4);
sig_vec=linspace(0.001,20,100)
plot(sig_vec,1./(2./sig_vec.*tanh(sig_vec/2)),'--r','LineWidth',4)

xlabel('\sigma');
ylabel('a_{\infty}');
legend(['Numerics       ';...
    'Sample specific';...
    'Statistical    '],...
    'Location', 'SouthEast');
axis ([0 20 0.9 M])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/a_inf.eps');

% %%
% figure;
% axes('FontSize',24);
% grid on
% hold on;
% plot(SMF_vec/M/2,-real(mobility));
% legend(num2str(sigma_vec'));
% xlabel('s/2M');
% ylabel('\mu')
% % axis([0 SMF_vec(end)/2/M, 0 .5])
%
%
%
% figure;
% axes('FontSize',24);
% grid on
% hold on;
% plot(SMF_vec/M/2,real(Diffusion));
% legend(num2str(sigma_vec'));
% xlabel('s/2M');
% ylabel('D')
% % axis([0 SMF_vec(end)/2/M, 0 .5])
