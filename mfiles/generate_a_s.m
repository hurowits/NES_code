clear all
DataPath='/Users/danielhurowitz/PROJ/NES/Figs/';
cols = 'rgbckym';c=0;
markers='*osd.'
M = 10;




s0 = (linspace(0,1,M)*2-1);
% s0 = randn(1,M);
for iR=1:1
    randInd1 = randperm(M);
    randInd2 = randperm(M);
    %     s0=rand(1,M)*2-1;
    % sigma=0.2;
    
    %         sigma_vec=[linspace(0.1,20,15)];
    %     sigma_vec=0.1
    %         sigma=sigma_vec(iS)   ;
    sigma_vec=3.5;
    %     SMF_vec = [linspace(0,M*3*sigma,300)];
    
    %         SMF_vec=1e3
    for iS=1:length(sigma_vec)
        sigma = sigma_vec(iS);
        %                 SMF_vec = logspace(-2,log10(M*3*sigma),1e3);
        SMF_vec=linspace(2/M,3*sigma*M,5*1e2);
        w2=zeros(1,length(SMF_vec));
        velocity=zeros(1,length(SMF_vec));
        Diffusion=zeros(1,length(SMF_vec));
        
        v_drda=zeros(1,length(SMF_vec));
        D_drda=zeros(1,length(SMF_vec));
        v_a=zeros(1,length(SMF_vec));
        D_a=zeros(1,length(SMF_vec));
        
        mu1=zeros(1,length(SMF_vec));
        mu2=zeros(1,length(SMF_vec));
        G=ones(1,M);
        for q=1:length(SMF_vec)
            
            SMF = SMF_vec(q);
            s=SMF/M;
            s_p = sigma*s0(randInd1)+SMF/M*ones(1,M);
            
            w_p = G.*exp(s_p/2);
            w_m = G.*exp(-s_p/2);
            
            w_harm1(q) = 1./mean(1./w_p);
            w_harm2(q) = 1./mean(1./w_p.^2);
            
            
            [velocity(q),Diffusion(q)] = VD(w_p,w_m);
            %             [velocity(iS),Diffusion(iS)] = VD(w_p,w_m);
%             if(Diffusion(q)>0)
            if(q==1)
                a_s_0=M;
            else
                a_s_0=a_s(q-1);
            end
            
            f = @(x)myfun(x,s,velocity(q)./Diffusion(q)/M);
            options=optimoptions('fsolve','Display','off','TolFun',1e-7,'TolX',1e-7);
            [a_s(q) fval(q),exitflag(q)]=fsolve(f,a_s_0,options);
            w2(q)=(1-mean((w_m./w_p).^2));
%             end
        end
       
    end
end
    %%
    %
    % figure;
    % axes('FontSize',24);
    % grid on
    % hold on;
    % plot(nu_vec, (real(velocity)),nu_vec,-SMF/M^2,'LineWidth',4);
    % xlabel('\nu [scaled]');
    % ylabel('V');
    % axtype(1)
    % axis([1e-7 5 0 10])
    % print(gcf, '-depsc2', '/users/physics/hurowits/MyWebSite/PROJ/NES/Figs/DvsS_500.eps');
    
%     
%     figure(10);
%     axes('FontSize',24);
%     grid on
%     hold on
%     
%     
%     hold on
%     plot(a_s.*SMF_vec/M,velocity./Diffusion.*a_s/2/M,[markers(iR),cols(iR)])
%     %     plot(SMF_vec/M,tanh(a_s.*SMF_vec/2/M),'--r','LineWidth',2)
%     xlabel('a_s\times s')
%     ylabel('a_s \times  f(s)')
%     
%     plot(SMF_vec/M,tanh(SMF_vec/2/M),'k','LineWidth',2)
%     % plot(SMF_vec/M,tanh(a_inf(end)*SMF_vec/2/M),'k','LineWidth',4)
%     axis([SMF_vec(1)/M SMF_vec(end)/M -0.01 1.2]);
%     %     print(gcf, '-depsc2', [DataPath,'vd_scaled.eps']);
%     
    
  


%%

