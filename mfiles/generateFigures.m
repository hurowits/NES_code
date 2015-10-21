nBins=100;
colors='rgbckym'
iE = 1;
Estrength_vec =  logspace(-5,3,1000);
sigma_vec = linspace(0,10,50);


for iS = [10,20,30,40,50]
    figure;
    axes('FontSize',24);
    grid on
    hold on;
    c=1;
    for iE = 1:200:1000
        sigmaI = std(CurrentEnsemble(:,iE,iS));
        mu = mean(CurrentEnsemble(:,iE,iS));
        I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
        
        [P,x] = hist(I,nBins);
        P = P/sum(P(:));
        
        P_filtered = conv(P,ones(1,10),'same')/10;
        
        %         plot(x,P_filtered,['.',colors(c)],'MarkerSize',8);
        plot(x,exp(-(x.^2)/2)/sum(exp(-(x.^2)/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    
    pleg=legend([repmat('\epsilon^2 = ',[5 1]),num2str((Estrength_vec(1:200:1000).^2'),'% 8.2e')]);
    %     legend('boxoff')
    set(pleg,'FontSize',16)
    
    c=1;
    for iE = 1:200:1000
        sigmaI = std(CurrentEnsemble(:,iE,iS));
        mu = mean(CurrentEnsemble(:,iE,iS));
        I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
        
        [P,x] = hist(I,nBins);
        P = P/sum(P(:));
        
        P_filtered = conv(P,ones(1,10),'same')/10;
        plot(x,P_filtered,['.',colors(c)],'MarkerSize',12);
        %         plot(x,exp(-x.^2/2)/sum(exp(-x.^2/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    xlabel(['I, \sigma = ',num2str(sigma_vec(iS))]);
    ylabel('P(I)');
    %     title(['\sigma = ',num2str(sigma_vec(iS))]);
    axis([-7, 7, 0, 0.06])
    print(gcf, '-depsc2', ['P_sig_',num2str((iS))])
    
end

%%
for iE = 1:200:1000
    figure;
    axes('FontSize',24);
    hold on;
    grid on;
    
    c=1;
    for iS = [10,20,30,40,50]
        sigmaI = std(CurrentEnsemble(:,iE,iS));
        mu = mean(CurrentEnsemble(:,iE,iS));
        I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
        
        [P,x] = hist(I,nBins);
        P = P/sum(P(:));
        
        P_filtered = conv(P,ones(1,10),'same')/10;
        
        %         plot(x,P_filtered,['.',colors(c)],'MarkerSize',8);
        plot(x,exp(-x.^2/2)/sum(exp(-x.^2/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    
    pleg= legend([repmat('\sigma = ',[5 1]),num2str((sigma_vec(  [10,20,30,40,50])'),'% 8.2e')]);
    %     legend('boxoff')
    set(pleg,'FontSize',16)
    
    c=1;
    for iS = [10,20,30,40,50]
        sigmaI = std(CurrentEnsemble(:,iE,iS));
        mu = mean(CurrentEnsemble(:,iE,iS));
        I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
        
        [P,x] = hist(I,nBins);
        P = P/sum(P(:));
        
        P_filtered = conv(P,ones(1,10),'same')/10;
        
        plot(x,P_filtered,['.',colors(c)],'MarkerSize',10);
        %         plot(x,exp(-x.^2/2)/sum(exp(-x.^2/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    %     title(['\epsilon^2 = ',num2str(Estrength_vec(iE)^2)]);
    axis([-7, 7, 0, 0.06])
    xlabel(['I, \epsilon^2 = ',num2str(Estrength_vec(iE)^2,'%8.2e')]);
    ylabel('P(I)');
    print(gcf, '-depsc2', ['P_eps_',num2str((iE))])
end
%%
sig = linspace(0,10,50);
% gmax = 2*exp(sig/2).*sinh(sig/2)./sig;
% gmin = 2*exp(-sig/2).*sinh(sig/2)./sig;
%  

gmax = 10.^(sig/2)./sig./log(10) .*(10.^(sig/2)-10.^(-sig/2));
gmin = 10.^(-sig/2)./sig./log(10) .*(10.^(sig/2)-10.^(-sig/2));
figure;
axes('FontSize',24);
loglog(sig,1./gmin,'r','LineWidth',2);
hold on;
grid on

plot(sig,1./gmax,'b','LineWidth',2);
for iE = 1:200:1000
    for iS = [10,20,30,40,50]
        
            if(iS==30&&iE==201)
                plot(sigma_vec(iS),Estrength_vec(iE).^2,'dr','MarkerSize',15);
            elseif(iS==30&&iE==601)
                plot(sigma_vec(iS),Estrength_vec(iE).^2,'sg','MarkerSize',15);
            elseif (iS==30&&iE==801)
                plot(sigma_vec(iS),Estrength_vec(iE).^2,'*b','MarkerSize',15);
            else
                plot(sigma_vec(iS),Estrength_vec(iE).^2,'ok','MarkerSize',8);
            end
            
    end
end
% set(gca,'XTick',log10(sigma-vec))
% set(gca,'YTick',log10(Estrength_vec.^2))
legend(['g_{min}^{-1}';'g_{max}^{-1}'],'Location','SouthWest');
xlabel('\sigma');
ylabel('\epsilon^2');
% print(gcf, '-depsc2', 'regions');
% axis([sig(1),sig(end),Estrength_vec(1)^2/10,Estrength_vec(end)^2]);
% print(gcf, '-depsc2', ['eps_sig_regions'])

