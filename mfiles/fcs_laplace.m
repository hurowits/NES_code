clear all;
% clear v D a_inf h_inf
load randInd.mat
% close all;
options = optimoptions('fsolve','Display','off');
flag=1;% flag=1 for g(lambda), flag=0 for (sigma,delta) dependence
numRealizations=1;
M=30;
% load randomIndices.mat
Delta_vec=linspace(0,5,100);[0,2,3];
sigma_vec=linspace(0,5,100);;[0,3,4];

% Delta_vec=[0,3,5];
% sigma_vec=[0,3,5];
% Delta_vec=5;3;[2]
% sigma_vec=0;3;
sigma_vec=0;3;0.5;4;
Delta_vec=0;3;0;
% sigma_vec=[0,3,4];
% Delta_vec=[0,3,4]
s0 = linspace(0,1,M)*2-1;
G0 = linspace(0,1,M)*2-1;
% SMF_vec=linspace(0,6,1e3)*M;M*[0.1,0.5,1,2,4];M/10;4*M;M/100;M;M/100;;M*4;
SMF_vec=[0.1*M,M,3*M,6*M,10*M];
SMF_vec=linspace(0,10*M,1e3);
% k_vec=linspace(0,2*pi,2*1024);

SMF_vec=[0,4*M];
k_vec = linspace(-130,130,1e3);
k_vec2 = linspace(0,2*pi,1e3);
% k_vec = k_vec(1:end-1);
colors='rgbckym'
if(flag==1)
    figure(1);
    axes('FontSize',24);
    grid on
    hold on
end

for iR=1:numRealizations
%     randInd1 = randperm(M);
%     randInd2 = randperm(M);
    
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
                        
                        [lambda_0(iK) ] = charFn( W, -1i*k);
                        %                     [lambda_02(iK) ] = charFn( W, -1i*(SMF-k));
                        %
                        %                     lambda_AA(iK) = (w_m(1)*exp(1i*k/M)+w_p(1)*exp(-1i*k/M))-(w_p(1)+w_m(1));
                        %                     lambda_AA_2(iK) = (w_m(1)*exp(-1i*k2)+w_p(1)*exp(1i*k2))-(w_p(1)+w_m(1));
                    end
                end
                h(iSMF,iSigma,iDelta)=charFn(W,-1i*SMF/2);
                
                delta_k=1e-2*(-1i) ;
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
                
%                 x_1(iR) = -v(iSMF,iSigma,iDelta)/(s*M*D(iSMF,iSigma,iDelta));
%                 x_2(iR) = real(h(iSMF,iSigma,iDelta))/(s*M*v(iSMF,iSigma,iDelta));
%                 
%                 S_1(iR,iSigma,iDelta) = abs(fsolve(@(x)2./x.*tanh(x/2)+v(iSMF,iSigma,iDelta)/(s*M*D(iSMF,iSigma,iDelta)),s*M,options));
%                 S_2(iR,iSigma,iDelta) = abs(fsolve(@(x)1./x.*(1-cosh(x/2))./sinh(x/2)+real(h(iSMF,iSigma,iDelta))/(s*M*v(iSMF,iSigma,iDelta)),s*M,options));
%                 
                %
                if(flag==1)
                    plot(k_vec/M,-real(lambda_0),colors(iSMF),'LineWidth',2);
                    plot(k_vec/M, -v(iSMF)*k_vec-D(iSMF)*k_vec.^2,['--',colors(iSMF)],'LineWidth',2);
                    xlabel('\lambda/M');
                    ylabel('g(\lambda)');
                    %                 textLegend=[textLegend;'Delta = ', num2str(Delta_vec(iDelta))]
                    %                             legend(['\Delta =',num2str(Delta_vec(1));...
                    %                                     'v\lambda+D\lambda^2'])
                    %                 '\Delta=',num2str(Delta_vec(2))]);
                    %
                    %                 title(['\sigma=',num2str(sigma),',s=',num2str(s)]);
                end
                
                
                
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
end
% legend(num2str(1/M*SMF_vec'));
axis([-4 4 -20 10])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/g_lambda_M.eps'])
%
%%
figure;
axes('FontSize',24);
grid on;
hold on;
plot(S_1/M,S_2/M,'.');
xlabel('S_1(v/sD)');
ylabel('S_2(h/sv)');
title(['\Delta=',num2str(Delta),', \sigma=',num2str(sigma),', s=',num2str(SMF/M)])
%%
% figure;plot(k_vec,(lambda_0),k_vec,lambda_02,'--')