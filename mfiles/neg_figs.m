for iSMF=1:length(SMF_vec)
    s=SMF_vec(iSMF)/M;
    %       figure(1);
    figure;
    axes('FontSize',24);
    grid on
    hold on
    imagesc(Delta_vec,sigma_vec,real(squeeze(4*h(iSMF,:,:)./(s*M*v(iSMF,:,:)))),[0 1]);axis image
   colorbar
    xlabel('\Delta');
    ylabel('\sigma');
    title(['4h/vs, s=',num2str(s),', M=',num2str(M)]);
    hcb=colorbar;
    set(hcb,'FontSize',24);
%   
%     print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/vhs_',num2str(iSMF),'.eps'])
    
end
%%
for iSMF=1:length(SMF_vec)
    s=SMF_vec(iSMF)/M;
    %       figure(1);
    figure;
    axes('FontSize',24);
    hold on
    %     axes('FontSize',24);%hold on;
%     contour(Delta_vec,sigma_vec,squeeze(-v(iSMF,:,:)./(D(iSMF,:,:)*s*M)),50);axis image
    imagesc(Delta_vec,sigma_vec,squeeze(-v(iSMF,:,:)./(D(iSMF,:,:)*s*M)),[0 1]);axis image
%     imagesc(Delta_vec,sigma_vec,squeeze(h(iSMF,:,:)));axis image

    %     imagesc(Delta_vec,sigma_vec,-1./squeeze(h(iSMF,:,:)./(v(iSMF,:,:).^2./(4*D(iSMF,:,:)))),[0.5 1.2]);axis image
    %             imagesc(Delta_vec,sigma_vec,squeeze(abs(h(iSMF,:,:))));axis image
    %     imagesc(Delta_vec,sigma_vec,squeeze(-v(iSMF,:,:)./(D(iSMF,:,:))));axis image;colorbar;
    colorbar
    xlabel('\Delta');
    ylabel('\sigma');
    title(['v/(Ds), s=',num2str(s),', M=',num2str(M)]);
    hcb=colorbar;
    set(hcb,'FontSize',24);
%     title(['D/(h/s), s=',num2str(s),', M=',num2str(M)]);
    %     title(['(v^2/4D)/h, s=',num2str(s),', M=',num2str(M)]);
    %         title(['h, s=',num2str(s),', M=',num2str(M)]);
    %     title(['v/D, s=',num2str(s),', M=',num2str(M)]);
%         print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/vDs_',num2str(iSMF),'.eps'])
    %     print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/hvD_s_',num2str(iSMF),'.eps'])
    %     print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/h_s_',num2str(iSMF),'.eps'])
    %     print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/vD_s_',num2str(iSMF),'.eps'])
    
end
%%
a = squeeze(a_inf(end,:,:));
for iSMF=1:length(SMF_vec)
    s=SMF_vec(iSMF)/M;
    %       figure(1);
    figure;
    axes('FontSize',24);
    hold on
    %     axes('FontSize',24);%hold on;
    imagesc(Delta_vec,sigma_vec,real(squeeze(-v(iSMF,:,:)./(D(iSMF,:,:)))./(2*tanh(s*a/2)./(a/M))),[0 1]);axis image
%     contour(Delta_vec.^2,sigma_vec.^2,a,100);axis image
    colorbar
    xlabel('\Delta');
    ylabel('\sigma');
    title([' s=',num2str(s),', M=',num2str(M)]);
    hcb=colorbar;
    set(hcb,'FontSize',24);
%   
%     title(['D/(h/s), s=',num2str(s),', M=',num2str(M)]);
    %     title(['(v^2/4D)/h, s=',num2str(s),', M=',num2str(M)]);
    %         title(['h, s=',num2str(s),', M=',num2str(M)]);
    %     title(['v/D, s=',num2str(s),', M=',num2str(M)]);
%         print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/vD_vDinf_',num2str(iSMF),'.eps'])
    
end
%%
a = squeeze(a_inf(end,:,:));
for iSMF=1:length(SMF_vec)
    s=SMF_vec(iSMF)/M;
    %       figure(1);
    figure;
    axes('FontSize',24);
    hold on
    %     axes('FontSize',24);%hold on;
    imagesc(Delta_vec,sigma_vec,-real(squeeze(h(iSMF,:,:)./h_inf(1,:,:))));axis image
%     imagesc(Delta_vec,sigma_vec,-real(squeeze(h(iSMF,:,:))));axis image
    colorbar
    xlabel('\Delta');
    ylabel('\sigma');
    title(['h/ h_{\infty} , s=',num2str(s),', M=',num2str(M)]);
%     title(['D/(h/s), s=',num2str(s),', M=',num2str(M)]);
    %     title(['(v^2/4D)/h, s=',num2str(s),', M=',num2str(M)]);
    %         title(['h, s=',num2str(s),', M=',num2str(M)]);
    %     title(['v/D, s=',num2str(s),', M=',num2str(M)]);
%         print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/h_hinf_',num2str(iSMF),'.eps'])
    
end
%%
x=SMF_vec/M;
figure;
axes('FontSize',24);
set(gca,'LineStyleOrder','-|-.|:')
hold on
grid on
txtLegend=[];
    for iS=1:length(sigma_vec)
            sigma=sigma_vec(iS);
        for iD=1:length(Delta_vec)
            Delta=Delta_vec(iD);
%             a_inf_ana(iS,iD) = sigma*Delta/4 *coth(sigma/2)*coth(Delta/2);
%         plot(x*a_inf(end,iS,iD),squeeze(-v(:,iS,iD)./D(:,iS,iD))/M./(x'),colors(iS),'LineWidth',4);
plot(a_inf(end,iS,iD)*x,a_inf(end,iS,iD)*squeeze(-v(:,iS,iD)./D(:,iS,iD))/M,colors(iS),'LineWidth',4);
%         plot(x*a_inf_ana(iS,iD),a_inf_ana(iS,iD)*squeeze(-v(:,iS,iD)./D(:,iS,iD))/M,'r','LineWidth',1);
        hold all
        txtLegend=   [txtLegend ;...
                      '\sigma=',num2str(sigma_vec(iS)),', \Delta=',num2str(Delta_vec(iD))]
    end
end
x=linspace(0,50,1e3);
plot(x,2*tanh(x/2),'--k','LineWidth',4);
% plot(x,2*tanh(M*x/2)/M,'k','LineWidth',4);

xlabel('sa_{\infty}');
ylabel('(v/D)a_{\infty}');
l=legend([txtLegend;...
    '\sigma=0, \Delta=0']);
set(l,'FontSize',14)
axis([0 27 0 2.2])
% legend(['\Delta=',num2str(Delta_vec(1)),', \sigma=',num2str(sigma_vec(1));...
%         '\Delta=',num2str(Delta_vec(2)),', \sigma=',num2str(sigma_vec(1));...
%         '\Delta=',num2str(Delta_vec(1)),', \sigma=',num2str(sigma_vec(2));...
%         '\Delta=',num2str(Delta_vec(2)),', \sigma=',num2str(sigma_vec(2))]);
%         print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/vD.eps'])

%%
figure;
axes('FontSize',24);
grid on
hold on
plot(k_vec,real(lambda_0),colors(iSMF),'LineWidth',2);
plot(k_vec, v(iSMF)*k_vec+D(iSMF)*k_vec.^2,['--',colors(iSMF)],'LineWidth',2);
xlabel('\lambda');
ylabel('g(\lambda)');
%