clear all
DataPath='/Users/danielhurowitz/PROJ/NES/Figs/';
cols = 'rgbckym';c=0;
markers='*osd.'
M = 50;


figure(11);

axes('FontSize',24);
grid on
hold on

s0 = (linspace(0,1,M)*2-1);

% s0 = randn(1,M);
randInd1 = randperm(M);
    randInd2 = randperm(M);
    Delta=2;
for iR=1:2
    
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
        %         SMF_vec = logspace(log10(1/M),log10(M*6),1e3);
        SMF_vec=linspace(0,8*M,1000);
        w2=zeros(1,length(SMF_vec));
        %         SMF_vec = [linspace(0,M*3*sigma,300)];
        
        %         SMF_vec = log10(M*2*sigma);
        %         SMF_vec=linspace(sigma/2*M,sigma*M,100);
        velocity=zeros(1,length(SMF_vec));
        Diffusion=zeros(1,length(SMF_vec));
        %         velocity=zeros(1,length(sigma_vec));
        %         Diffusion=zeros(1,length(sigma_vec));
        %
        
        SMF=zeros(1,length(SMF_vec));
        v_drda=zeros(1,length(SMF_vec));
        D_drda=zeros(1,length(SMF_vec));
        v_a=zeros(1,length(SMF_vec));
        D_a=zeros(1,length(SMF_vec));
        a_s=M*ones(1,length(SMF_vec));
           a_s2=M*ones(1,length(SMF_vec));
        mu1=zeros(1,length(SMF_vec));
        mu2=zeros(1,length(SMF_vec));
        if(iR==1)
            G=ones(1,M);
        else
            G = exp(rand(1,M)*Delta);
        end
        for q=1:length(SMF_vec)
            
            SMF = SMF_vec(q);
            s=SMF/M;
            s_p = sigma*s0(randInd1)+SMF/M*ones(1,M);
            
            w_p = G.*exp(s_p/2);
            w_m = G.*exp(-s_p/2);
            
            w_harm1(q) = 1./mean(1./w_p);
            w_harm2(q) = 1./mean(1./w_p.^2);
            a = 2/sigma*exp(-s/2)*sinh(sigma/2);
            b = exp(-s)/sigma * sinh(sigma);
            c = 1/2/sigma * exp(-2*s) * sinh(2*sigma);
            d = exp(-s)/sigma * sinh(sigma);
            e = 2/3/sigma * exp(-3/2*s) * sinh(3*sigma/2);
            
            
            [velocity(q),Diffusion(q)] = VD(w_p,w_m);
            %             [velocity(iS),Diffusion(iS)] = VD(w_p,w_m);
            %             if(Diffusion(q)>0)
            %                 if(q==1)
            %                     a_s_0=M;
            %                 else
            %                     a_s_0=a_s(q-1);
            %                 end
            %                 f = @(x)myfun(x,s,velocity(q)./(Diffusion(q))/M);
            %                 options=optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',10000,'MaxFunEvals',1000);
            %                 [a_s(q) fval(q) exitflag(q)]=fsolve(f,a_s_0,options);
            %                 if(a_s(q)>M)
            %                     a_s(q)=M;
            %                 end
            % %                 if (exitflag(q)==-2)
            % %                     break;
            % %                 end
            %             end
            
            w2(q)=(1-mean((w_m./w_p).^2));
            if (SMF>=0)
                mu1(q) = mean(w_m./w_p);
                mu2(q) = mean((w_m./w_p).^2);
                mu05(q) = mean((w_m./w_p).^0.5);
                if (mu1(q)<1)
                    v_drda(q) = (1- mean(w_m./w_p))./(mean(1./w_p));
                    V_a(q) =  (1-b)/a;
                else
                    
                    v_drda(q) = 0;
                    V_a(q) = 0;
                end
                
                if (b<1)
                    V_a(q) =  (1-b)/a;
                else
                    
                    V_a(q) = 0;
                end
                
                if(mu2(q)<1)
                    D_drda(q) = 0.5*(mean(1./w_p)).^(-3)*...
                        (1-(mean(w_m./w_p)).^2)./...
                        (1-mean((w_m./w_p).^2))*...
                        (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
                        2*mean(w_m./w_p.^2)*mean(1./w_p));
                    
                    %         end
                    %                     if(mu1(q)>1)
                    %                         D=0;
                    %                     end
                else
                    D_drda(q)=Inf;
                end
%                 if(mu05(q)>1)
%                     D_drda(q)=0;
%                 end
                if(c<1)
                    
                    D_a(q) = 1/2 * a^(-3)*(1 - b^2)/(1 - c)*(d*(1 - b) + 2*e*a);
                    
                    %         end
                    if(b>1)
                        D_a(q)=0;
                    end
                else
                    
                    D_a(q)=Inf;
                end
                %         mu1(q) = mean(w_m./w_p);
                %         mu2(q) = mean((w_m./w_p).^2);
                %
            else
                mu1(q) = mean(w_p./w_m);
                mu2(q) = mean((w_p./w_m).^2);
                if(mu1(q)<1)
                    v_drda(q) = -(1- mean(w_p./w_m))./(mean(1./w_m));
                else
                    v_drda(q)=0;
                end
                D_drda(q) = 0.5*(mean(1./w_m)).^(-3)*...
                    (1-(mean(w_p./w_m)).^2)./...
                    (1-mean((w_p./w_m).^2))*...
                    (mean(1./w_m.^2)*(1-mean(w_p./w_m))+...
                    2*mean(w_p./w_m.^2)*mean(1./w_m));
                %
                %         mu1(q) = mean(w_p./w_m);
                %         mu2(q) = mean((w_p./w_m).^2);
                
                %         v_drda(q,iSig) = (1- mean(w_m./w_p))./(mean(1./w_p));
                %         D_drda(q,iSig) = 0.5*(mean(1./w_p)).^(-3)*...
                %             (1-(mean(w_m./w_p))f.^2)./...
                %             (1-mean((w_m./w_p).^2))*...
                %             (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
                %             2*mean(w_m./w_p.^2)*mean(1./w_p));
                %
                
            end
            if (D_drda(q)<0)
                D_drda(q)=Inf;
            end
            
            
        end %for q
        a_inf(iS) = mean((1./w_p).^2)/mean(1./w_p)^2;
        a_inf1(iS) = 2*Diffusion(end)./velocity(end);
        vd_smooth = conv((velocity./Diffusion),ones(10,1),'same')/10;
        for q=1:length(SMF_vec)
            if(Diffusion(q)>0)
                if(q==1)
                    a_s_0=M;
                else
                    a_s_0=a_s2(q-1);
                end
                a_s_0=1;
                f = @(x)myfun(x,s,vd_smooth(q)/M);
                options=optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-9,'MaxIter',100000,'MaxFunEvals',10000);
                [a_s2(q) fval(q) exitflag(q)]=fsolve(f,a_s_0,options);
                if(a_s2(q)>M)
                    a_s2(q)=M;
                end
                %                 if (exitflag(q)==-2)
                %                     break;
                %                 end
            end
        end
        
    end
    
%%
figure(11);
% axes('FontSize',24);
% grid on
% hold on

plot(SMF_vec/M, SMF_vec/M/2,'--b','LineWidth',4);
%     plot(SMF_vec/M, (V_a./D_a)/2,'-.r','LineWidth',4);
% plot(SMF_vec/M, velocity./Diffusion/M/2,[markers(iR),cols(iR)]);
plot(SMF_vec/M, velocity./Diffusion/M/2,'r','LineWidth',2);
plot(SMF_vec/M,tanh(SMF_vec/M/2),'k','LineWidth',4)

plot(SMF_vec/M,1/a_inf(end)*tanh(a_inf(end)*SMF_vec/M/2),'k','LineWidth',4)
plot(SMF_vec/M,1/M*tanh(M*SMF_vec/M/2),'k','LineWidth',4)
plot(SMF_vec/M, (v_drda./D_drda)/2,'--b','LineWidth',4);
% plot(SMF_vec/M, 1./a_s.*tanh(a_s.*SMF_vec/M/2),'--g','LineWidth',4);
%     plot(SMF_vec/M, (velocity./Diffusion)/M/2.*tanh(SMF_vec/2),'--b','LineWidth',4);
axis([SMF_vec(1)/M SMF_vec(end)/M -0.01 1]);
% 
% plot([1/M 1/M], [1e-5 40],'--k')
% plot([s05 s05],[1e-5 40],'--k')
% plot([s1 s1],[1e-5 40],'--k')
% plot([s2 s2],[1e-5 40],'--k')
% 
xlabel('s');
ylabel('v/2D');
% hold off

%         legend(['ESR       ';...
%
%         %         'N \rightarrow \infty';...
%         'numerics  ';...
%         'a_0       ';...
%         'a_{\infty}';...
%         'Na_0      ';...
%         'Derrida   '],'Location','SouthEast');
%     print(gcf, '-depsc2', [DataPath,'v_D_rat.eps']);

    
end

%%

a_s_N=2*M*Diffusion./velocity.*tanh(SMF_vec/2);
figure;
axes('FontSize',24);
grid on
hold on
plot(SMF_vec/M,a_s,'r','LineWidth',4);
plot(SMF_vec/M,a_inf(end)./w2,'--b','LineWidth',4);
% plot(SMF_vec/M,a_s_N,'--k','LineWidth',4);


xlabel('s');
ylabel('a_s');
% axtype(3);
legend(['N=40    ';
    'N=\infty']);
axis([.009 SMF_vec(end)/M 0 1.1*M]);
% print(gcf, '-depsc2', [DataPath,'a_s.eps']);

%%
s=SMF_vec/M;
figure;
axes('FontSize',24);
hold on
plot(SMF_vec/M, D_a,'b','LineWidth',4)
plot(SMF_vec/M, V_a.*D_a(end)/V_a(end)./(1-exp(-2*s)*sinh(2*sigma)/2/sigma),'g','LineWidth',4)
legend(['D                       ';...
    'a_{inf} v /(1-<(w/w)^2>)']);
xlabel('s')
ylabel('D')
axis([0 6 1 25])
% print(gcf, '-depsc2', [DataPath,'D_approx.eps']);
% axis([0 12 0 8])
% plot(SMF_vec/M,velocity.*a_s*M,'r')
% plot(SMF_vec/M,Diffusion*M^2,'c')


%%
f = @(x)myfun3(x,sigma);
options=optimoptions('fsolve','Display','off','TolFun',1e-7,'TolX',1e-7);
[s05]=fsolve(f,1,options);

f = @(x)myfun1(x,sigma);
options=optimoptions('fsolve','Display','off','TolFun',1e-7,'TolX',1e-7);
[s1 ]=fsolve(f,1,options);
f = @(x)myfun2(x,sigma);
options=optimoptions('fsolve','Display','off','TolFun',1e-7,'TolX',1e-7);
[s2 ]=fsolve(f,1,options);

%
% sp1 = spmak([0 SMF_vec/M],mu1-1);
% s1 = fnzeros(sp1);
%
%
% sp2 = spmak([0 SMF_vec/M],mu2-1);
% s2 = fnzeros(sp2)

figure;
axes('FontSize',24);
% grid on
hold on

plot(SMF_vec(1:end)/M,velocity(1:end)*M./w_harm1(1:end),'r','LineWidth',4);
plot(SMF_vec/M,Diffusion*M^2./w_harm1,'b','LineWidth',4);
plot(SMF_vec(SMF_vec>1)/M,a_s((SMF_vec>1))/M,'g','LineWidth',4);
plot(SMF_300/300,a_s_300/300,'--g','LineWidth',3);
plot(SMF_300/300,300*v_300./w_harm1_300,'--r','LineWidth',3)

% plot(SMF_vec/M,D_drda./w_harm1,':b','LineWidth',4)
% plot(SMF_300/300,conv(a_s_300/300,ones(20,1),'same')/20,'-g','LineWidth',2)
% plot(SMF_vec/M,a_s.*velocity./(w_harm1)*M/2)

% plot([1/M 1/M], [0 4],'--k')
plot([s05 s05],[0 4],'--k')
plot([s1 s1],[0 4],'--k')
plot([s2 s2],[0 4],'--k')


% plot(a_s_drda.*SMF_vec/M,y1d,'--r','LineWidth',4);
% plot(a_s_drda.*SMF_vec/M,y2d,'--b','LineWidth',4);
%

xlabel('s');

% hleg1=legend(['v(s)  ';...
%     'D(s)  ';...
%     'a(s)/N'],'Location','SouthEast');
% set(hleg1,'Interpreter','tex');
% axis([1e-2 SMF_vec(end)/M 1e-8 1e-3]);
axis([0 8  0 1.6]);
% axtype(1);
text(s05-0.7,1.9,'s_{1/2}','FontSize',22)
text(s1-0.4,1.9,'s_1','FontSize',22)
text(s2-0.4,1.9,'s_2','FontSize',22)

text(4,1.3, 'D(s)','FontSize',18)
text(3,0.6, 'v(s)','FontSize',18)

ax2=axes('FontSize',18);

% ax2=axes;
hold on
% grid on
plot(SMF_vec(2:end-1)/M,velocity(2:end-1)*M./w_harm1(2:end-1),'r','LineWidth',4);
plot(SMF_vec/M,Diffusion*M^2./w_harm1,'b','LineWidth',4);

plot([1/M 1/M], [1e-7 1e0],'--k')
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
text(1.1/M,3*10^(-1),'1/N','FontSize',14)

axis([1/M/10 1 3*1e-7 1e0])
axtype(3)
set(ax2,'units','normalized','position',[.62 .25 .25 .25])
set(ax2,'box','on')

xlabel('s');
set(get(gca,'XLabel'),'Position',[.8 3*10^(-5.8)])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
% print(gcf, '-depsc2', [DataPath,'vda.eps']);

% axtype(3);
% axis([SMF_vec(1)/M SMF_vec(end)/M -.1 a_s(1)]);
%%


%%
figure;
% hold on
axes('FontSize',24);
grid on;
hold on;
%     plot(sigma_vec,M./velocity./Diffusion/2,'ro','MarkerSize',10,'LineWidth',4);
plot(sigma_vec,a_inf,'ro','MarkerSize',10,'LineWidth',4);
% plot(sigma_vec,a_inf,'rs','MarkerSize',10,'LineWidth',4);
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

sig_vec=linspace(0.001,20,100);
plot(sig_vec,1./(2./sig_vec.*tanh(sig_vec/2)),'--r','LineWidth',4)
plot(sig_vec,1./(log(sinh(sig_vec)./sig_vec)),'k:','LineWidth',2);
xlabel('\sigma');
ylabel('a_{\infty}');
legend(['Numerics       ';...
    'Sample specific';...
    'Statistical    '],...
    'Location', 'SouthEast');
axis ([0 20 0 M])
%%
%
% figure;
% axes('FontSize',24);
% grid on
% hold on;
% plot(SMF_vec/M,abs(velocity)*M,'-g','LineWidth',4);
%
% plot(SMF_vec/M,v_drda,'-b','LineWidth',4);
% plot(SMF_vec/M,V_a,'--r','LineWidth',4)
% % plot(nu_vec,SMF/M,'r','LineWidth',4);
% % axtype(1)
% xlabel('s')
% ylabel('v');
% legend(['numerical';...
%         'derrida  '],'Location','NorthWest');
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/brm_figs/v_derrida.eps')
%
%     figure;
%     axes('FontSize',24);
%     grid on
%     hold on;
%     plot(SMF_vec/M,D_drda,'b','LineWidth',4)
%     plot(SMF_vec/M,abs(Diffusion)*M^2,'g','LineWidth',4);
%     plot(SMF_vec/M,D_a,'--r','LineWidth',4);
%     % plot(nu_vec,abs(SMF)/M,'r','LineWidth',4);
%     xlabel('s');
%     ylabel('D')
%     % axtype(3);
%     legend(['derrida ';...
%         'numerics']);
%     title(['\sigma=',num2str(sigma),' N=',num2str(M)]);
%     axis([SMF_vec(1)/M SMF_vec(end)/M 0 10])
% % %     print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/brm_figs/D_derrida.eps')



%%
% figure;
% axes('FontSize',24);
% grid on
% hold on;
% plot(SMF_vec/M,log(mu1),'--b','LineWidth',4);
% plot(SMF_vec/M,log(mu2),'g','LineWidth',4);
% % plot(nu_vec,SMF/M,'r','LineWidth',4)
% % plot(SMF_vec/M,ones(size(nu_vec)),'--k');
% xlabel('\nu')
% legend(['\mu_1';...
%     '\mu_2']);
% axtype(1);
% axis([nu_vec(1) nu_vec(end) -.5 1.9]);
%     print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/brm_figs/derrida_mu.eps')
%%
%
% figure;
% grid on
% hold on;
% plot(nu_vec,v_drda,'b','LineWidth',4)
% plot(nu_vec,-velocity*M,'--b','LineWidth',4)
% plot(nu_vec,(D_drda/M^2),'r','LineWidth',4)
% plot(nu_vec,(Diffusion),'--r','LineWidth',4)
% plot(nu_vec,mu1,'-.b','LineWidth',4);
% plot(nu_vec,mu2,'-.r','LineWidth',4);
% axtype(1)