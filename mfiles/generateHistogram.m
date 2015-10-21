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

%% Estimate P(I;epsilon^2,sigma)

nBins=100;
colors='rgbckym'
iE = 1;
Estrength_vec =  logspace(-5,3,1000);
sigma_vec = linspace(0,10,50);
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
    legend([repmat('\sigma = ',[5 1]),num2str((sigma_vec(  [10,20,30,40,50])'),'% 8.2e')]);

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
    title(['\epsilon^2 = ',num2str(Estrength_vec(iE)^2)]);
    print(gcf, '-depsc2', ['P_eps_',num2str(Estrength_vec(iE).^2)])
end
%%
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
        plot(x,exp(-x.^2/2)/sum(exp(-x.^2/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    legend([repmat('\epsilon^2 = ',[5 1]),num2str((Estrength_vec(1:200:1000).^2'),'% 8.2e')]);
    c=1;
    for iE = 1:200:1000
        sigmaI = std(CurrentEnsemble(:,iE,iS));
        mu = mean(CurrentEnsemble(:,iE,iS));
        I = (CurrentEnsemble(:,iE,iS) - mu)/sigmaI;
        
        [P,x] = hist(I,nBins);
        P = P/sum(P(:));
        
        P_filtered = conv(P,ones(1,10),'same')/10;
        plot(x,P_filtered,['.',colors(c)],'MarkerSize',8);
%         plot(x,exp(-x.^2/2)/sum(exp(-x.^2/2)),['-',colors(c)],'LineWidth',2);
        c=c+1;
        
    end
    title(['\sigma = ',num2str(sigma_vec(iS))]);
    print(gcf, '-depsc2', ['P_sig_',num2str(sigma_vec(iS))])

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
