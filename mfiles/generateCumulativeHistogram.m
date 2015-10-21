%% Load dataset
% pathDatasets = '/users/physics/hurowits/MyWebSite/PROJ/NES/NES_code/datasets/';
pathDatasets = '/Users/danielhurowitz/Desktop/NES/NES_code/datasets/';
compList = ls(pathDatasets);
compList=reshape(compList',[7,length(compList)/7])';
CurrentEnsemble = zeros(8054,1000,50);
iR = 1;
for iComp = 1:length(compList)
    currentPath = [pathDatasets,compList(iComp,1:6),'/'];
    fileList = dir(currentPath);
    for iFile = 3:length(fileList);
        I = load([currentPath,fileList(iFile).name],'Current');
        CurrentEnsemble(iR,:,:) = I.Current;
        iR = iR+1;
    end
    
end
load([currentPath,fileList(iFile).name]);


%%


%% Estimate P(I;epsilon^2,sigma)

nBins=100;
colors='rgbckym'
iE = 1;
Estrength_vec =  logspace(-5,3,1000);
sigma_vec = linspace(0,10,50);
for iS =30; [201,601,801]
    figure(1);
    axes('FontSize',24);
    hold on;
    grid on;
    figure(2);
    axes('FontSize',24);
    hold on;
    grid on;
    
    c=1;
    
    
    for iE = [201,601,801]
        
%         sigmaI = std(CurrentEnsemble(:,iE,iS));
%         mu = mean(CurrentEnsemble(:,iE,iS));
%         I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
%         
        I = log(abs(CurrentEnsemble(:,iE,iS)));
        I = (I- mean(I))/std(I);
          x = linspace(min(I),max(I),length(I));
        figure(1);
        plot(x,log(normcdf(x,0,1)),['-',colors(c)],'LineWidth',2);
        figure(2);
        plot(x,log(1-normcdf(x,0,1)),['-',colors(c)],'LineWidth',2);
        

        c=c+1;
        
    end
    c=1;
      figure(1);
    title(['\sigma = ',num2str(sigma_vec(iS))]);
    legend([repmat('\epsilon^2 = ',[3 1]),num2str((Estrength_vec( [201,601,801]).^2'),'% 8.2e')],'Location','SouthEast');
    xlabel('(log(|I|)-<log(|I|)>)/std(log(|I|))');
    ylabel('log(CDF)');
    figure(2);
    title(['\sigma = ',num2str(sigma_vec(iS))]);
    legend([repmat('\epsilon^2 = ',[3 1]),num2str((Estrength_vec( [201,601,801]).^2'),'% 8.2e')],'Location','SouthWest');
    xlabel('(log(|I|)-<log(|I|)>)/std(log(|I|))');
    ylabel('log(CDF^{-1})');
    for iE = [201,601,801]
        
%         sigmaI = std(CurrentEnsemble(:,iE,iS));
%         mu = mean(CurrentEnsemble(:,iE,iS));
%         I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
%         
        I = log(abs(CurrentEnsemble(:,iE,iS)));
        I = (I- mean(I))/std(I);
        figure(1);
        plot(sort(I),log((1:length(I))/length(I)),['o',colors(c)],'MarkerSize',5);
        figure(2);
        plot(sort(I),log((length(I)-(1:length(I)))/length(I)),['o',colors(c)],'MarkerSize',5);
        

        c=c+1;
        
    end
    
  
    %     print(gcf, '-depsc2', ['P_eps_',num2str(Estrength_vec(iE).^2)])
end
%%
nBins=100;
colors='rgbckym'
iE = 1;
Estrength_vec =  logspace(-5,3,1000);
sigma_vec = linspace(0,10,50);
for iS =20 % [10,20,30,40,50]
    
    figure;
    axes('FontSize',24);
    grid on
    hold on;
    c=1;
    for iE = 201 %1:200:1000
%         sigmaI = std(CurrentEnsemble(:,iE,iS));
%         mu = mean(CurrentEnsemble(:,iE,iS));
%         I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
%         
%         %         Ipos = I(I>0);
%         %         Ineg = I(I<0);
%         %         I = log(abs(I));
%         %         mu = mean(I); sigmaI=std(I);
%         %         plot((sort(I)),(normcdf(sort(I),mu,sigmaI)),['-',colors(c)],'LineWidth',2);
%         %         plot((sort(I)),((1-normcdf(sort(I),mu,sigmaI))),['-.',colors(c)],'LineWidth',2);
%         plot((sort(abs(I))),log10(normcdf(sort(I),0,1)),['-',colors(c)],'LineWidth',2);
%         plot((sort(abs(I))),((1-normcdf(sort(I),0,1))),['-.',colors(c)],'LineWidth',2);
%         
   
        I = log(abs(CurrentEnsemble(:,iE,iS)));
        I = (I- mean(I))/std(I);
        x = linspace(min(I),max(I),length(I));
        
%         plot(x,log10(normcdf(x,0,1)),['-','r'],'LineWidth',2);
        semilogy(x,1-normcdf(x,0,1),['-.','b'],'LineWidth',2);

        c=c+1;
        
    end
    %     legend([repmat('\epsilon^2 = ',[5 1]),num2str((Estrength_vec(1:200:1000).^2'),'% 8.2e')]);
    c=1;
    for iE = 201 % 1:200:1000
        sigmaI = std(CurrentEnsemble(:,iE,iS));
%         mu = mean(CurrentEnsemble(:,iE,iS));
%         I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
%         
        %          I = log(abs(I));
        %         mu = mean(I); sigmaI=std(I);
        %
%         %
%                plot((sort(abs(I))),log10((1:length(I))/length(I)),['.',colors(c)],'MarkerSize',10);
%                plot((sort(abs(I))),(length(I)-(1:length(I)))/length(I),['o',colors(c)],'MarkerSize',5);
%         
        %
        %        plot((sort(I)),(exp(-exp(-(sort(I)-mu)/sigmaI))),['-',colors(c+1)],'LineWidth',2);
        %        plot((sort(I)),1-(exp(-exp(-(sort(I)-mu)/sigmaI))),['-.',colors(c+1)],'LineWidth',2);
        %
        I = log(abs(CurrentEnsemble(:,iE,iS)));
        I = (I- mean(I))/std(I);
        
%         plot(sort(I),log10((1:length(I))/length(I)),['.','r'],'MarkerSize',10);
        semilogy(sort(I),(length(I)-(1:length(I)))/length(I),['o','b'],'MarkerSize',5);
%         
%         plot(linspace(min(I),max(I),length(I)),(1/sigma_vec(iS)*log(linspace(min(I),max(I),length(I)))))
%         
%         xlabel('log(|I|)');
%         ylabel('log10(CDF)');
% %        
        c=c+1;
        
    end
    legend(['log10(CDF)      ';...
            'CDF^{-1 }       ';...
            'log10(hist)     ';...
            'hist^{-1}       '])
    title(['\sigma = ',num2str(sigma_vec(iS)),' \epsilon^2 = ',num2str(Estrength_vec(iE)^2)]);
    %     print(gcf, '-depsc2', ['P_sig_',num2str(sigma_vec(iS))])
    
end

%% Estimate P(I;epsilon^2,sigma)
nBins=100;
P2 = zeros(nBins,1000,50);
x2 = zeros(nBins,1000,50);
P_est_image= zeros(nBins,1000,50);
for iE = 1:1000
    for iS = 5:50
        I = CurrentEnsemble(:,iE,iS);
        I = (I-mean(I))/std(I);
        [P2(:,iE,iS),x2(:,iE,iS)] = hist(I(:),nBins);
        P_est_image(:,iE,iS) = exp(-x2(:,iE,iS).^2/2)/sum(exp(-x2(:,iE,iS).^2/2));
        
    end
end
%%
for iS=[10,30,50]
    figure;
    axes('FontSize',24);
    imagesc(log10(Estrength_vec.^2),linspace(min(CurrentEnsemble(:)),max(CurrentEnsemble(:)),nBins),...
        (conv2(P2(:,:,iS),ones(10,1),'same')/10/8055));
    
    xlabel(['\epsilon^2, \sigma = ',num2str(sigma_vec(iS))]);
    ylabel('(I-<I>)/std(I)');
    colorbar;
    %  print(gcf, '-depsc2', ['P2_image_sig_',num2str(sigma_vec(iS))])
    figure;
    axes('FontSize',24);
    imagesc(log10(Estrength_vec.^2),linspace(min(CurrentEnsemble(:)),max(CurrentEnsemble(:)),nBins),...
        P_est_image(:,:,iS));
    
    xlabel(['\epsilon^2, \sigma = ',num2str(sigma_vec(iS))]);
    ylabel('(I-<I>)/std(I)');
    colorbar;
    %   print(gcf, '-depsc2', ['P_est_image_sig_',num2str(sigma_vec(iS))])
    
    figure;imagesc(log10(Estrength_vec.^2),linspace(min(CurrentEnsemble(:)),max(CurrentEnsemble(:)),nBins),...
        abs(P_est_image(:,:,iS)- (conv2(P2(:,:,iS),ones(10,1),'same')/10/8055)));
    colorbar;
    %        print(gcf, '-depsc2', ['dP_image_sig_',num2str(sigma_vec(iS))])
    
end
