clear all;
% close all;
flag=1;% flag=1 for g(lambda), flag=0 for (sigma,delta) dependence

M=30;
randInd1 = randperm(M);
randInd2 = randperm(M);
% load randomIndices.mat
Delta_vec=linspace(0,4,100);[0,3,5];
sigma_vec=linspace(0,4,100);;[0,3,5];

% Delta_vec=[0,3,5];
% sigma_vec=[0,3,5];
Delta_vec=0;3;[2]
sigma_vec=3;3;
% sigma_vec=0;3;0.5;4;
% Delta_vec=0;3;0;
s0 = linspace(0,1,M)*2-1;
G0 = linspace(0,1,M)*2-1;
SMF_vec=linspace(0,6,1e3)*M;M*[0.1,0.5,1,2,4];M/10;4*M;M/100;M;M/100;;M*4;
SMF_vec=M;M;[0.1*M,M,2*M,3*M,6*M,10*M];
%%


k_vec = (0:1023)*2*pi/1024; linspace(0,2*pi,1024);
% k_vec = k_vec(1:end-1);
colors='rgbckym'
if(flag==1)
            figure(1);
            axes('FontSize',24);
            grid on
            hold on
        end
for iSMF = 1:length(SMF_vec)
    SMF = SMF_vec(iSMF);
    textLegend=[];
    for iSigma = 1:length(sigma_vec)
        
        sigma=sigma_vec(iSigma);
        for iDelta=1:length(Delta_vec)
            Delta = Delta_vec(iDelta);
            
            G = exp(G0(randInd2)*Delta);
            
            
              s=SMF/M;
            s_p = sigma*s0(randInd1)+SMF/M*ones(1,M);
            
            w_p = G.*exp(s_p/2);
            w_m = G.*exp(-s_p/2);
            
            W = -diag(w_p + w_m([M,1:M-1]),0) + ...
                diag(w_m(1:M-1),1)+diag(w_p(1:M-1),-1);
            W(1,M) = w_p(M);
            W(M,1) = w_m(M);
            
            if(flag==1)
                for iK=1:length(k_vec)
                    k=k_vec(iK);
                    %                     k2=k_vec2(iK);
                    %                     [P_k(iK,:)] = charFn_fcs(W,-1i*k);
                    
                    [lambda_0(iK) ] = charFn( W, k);
                    [P(iK,:)] = charFn_fcs( W, k);
                        
%                     [lambda_02(iK) ] = charFn( W, -1i*(SMF-k));
                    %
                    %                     lambda_AA(iK) = (w_m(1)*exp(1i*k/M)+w_p(1)*exp(-1i*k/M))-(w_p(1)+w_m(1));
                    %                     lambda_AA_2(iK) = (w_m(1)*exp(-1i*k2)+w_p(1)*exp(1i*k2))-(w_p(1)+w_m(1));
                end
                P_f_bar = fftshift(mean(P,2));
                P_n =  fftshift(ifft(fftshift(P_f_bar)));
            end
            h(iSMF,iSigma,iDelta)=charFn(W,SMF/2);
            %min(lambda_0);
            
            delta_k=1e-2 ;
            lambda1_m = charFn(W,-delta_k/2);
            lambda1_p = charFn(W,delta_k/2);
            lambda2_m = charFn(W,-delta_k);
            lambda2_p = charFn(W,delta_k);
            lambda0 = 0;
            dlambdadk = -1i*(lambda1_p - lambda1_m)/delta_k;
            
            v(iSMF,iSigma,iDelta) = dlambdadk;
            D(iSMF,iSigma,iDelta) = -0.5*(lambda2_p + lambda2_m )/delta_k^2;
            a_inf(iSMF,iSigma,iDelta) = mean((1./w_p).^2)/mean(1./w_p)^2;
            h_inf(1,iSigma,iDelta) = min(w_p);
           
            %
            if(flag==1)
%                 plot(k_vec,(abs((lambda_0))),colors(iSMF),'LineWidth',2);
                plot(k_vec,((((lambda_0-1i*pi*k_vec*s/2+1i*pi/2)))),colors(iSMF),'LineWidth',2);
                
                %             plot(-1i*k_vec,real(lambda_AA),['--',colors(iDelta+1)]);
                %             plot(k_vec, -v*k_vec+D*k_vec.^2);
%                 plot(k_vec, v(iSMF)*k_vec+D(iSMF)*k_vec.^2,['--',colors(iSMF)],'LineWidth',2);
                %             plot(-1i*k_vec,real(lambda_AA_2),['--',colors(iDelta+2)]);
                %                 plot(k_vec,angle(lambda_0),['--',colors(iR)]);
                xlabel('\lambda/M');
                ylabel('g(\lambda)');
%                 textLegend=[textLegend;'Delta = ', num2str(Delta_vec(iDelta))]
%                             legend(['\Delta =',num2str(Delta_vec(1));...
%                                     'v\lambda+D\lambda^2'])
                %                 '\Delta=',num2str(Delta_vec(2))]);
                %
%                 title(['\sigma=',num2str(sigma),',s=',num2str(s)]);
            end
            
            
            %     P_x = ifft(exp(lambda_0));
            %     P_x = sum(fftshift(ifft((P_k)/2048)),2);
            %     P_x = P_x./sum(P_x);
            %     figure(3);
            %     hold on;
            %     plot(linspace(0,M,length(k_vec)),((real(P_x))),colors(iR));
            %     % figure(2);hold on;
            % plot(linspace(-M,M,length(k_vec)),angle((P_x)),colors(iR));
            %
        end
    end
%     %%
%     figure;
%     axes('FontSize',24);hold on;
%     imagesc(Delta_vec,sigma_vec,(-(v./D)));axis image
%     xlabel('\Delta');
%     ylabel('\sigma');
%     title(['v/D, s=',num2str(s),', M=',num2str(M)]);
%     %
%     figure;
%     axes('FontSize',24);hold on;
%     imagesc(Delta_vec,sigma_vec,(-h));axis image
%     xlabel('\Delta');
%     ylabel('\sigma');
%     title(['h, s=',num2str(s),', M=',num2str(M)]);
%     %
%     figure;
%     axes('FontSize',24);hold on;
%     imagesc(Delta_vec,sigma_vec,1./squeeze(-h(iSMF,:,:)./v(iSMF,:,:).^2./(4*D(iSMF,:,:))));axis image
%     xlabel('\Delta');
%     ylabel('\sigma');
%     title(['h/(v^2/4D), s=',num2str(s),', M=',num2str(M)]);
%     
end
% legend(num2str(1/M*SMF_vec'));

% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/g_lambda_M.eps'])
%
%%
% figure;plot(k_vec,(lambda_0),k_vec,lambda_02,'--')