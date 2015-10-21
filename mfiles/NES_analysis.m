%% Load dataset
% pathDatasets = '/users/physics/hurowits/MyWebSite/PROJ/NES/NES_code/datasets/';
pathDatasets = '/Users/danielhurowitz/Desktop/NES/NES_code/datasets/';
compList = ls(pathDatasets);
compList=reshape(compList',[7,length(compList)/7])';
CurrentEnsemble = zeros(8054,1000,50);
SMFEnsemble = zeros(8054,1000,50);
iR = 1;
for iComp = 1:length(compList)
    currentPath = [pathDatasets,compList(iComp,1:6),'/'];
    fileList = dir(currentPath);
    for iFile = 3:length(fileList);
        I = load([currentPath,fileList(iFile).name],'Current');
        SMF = load([currentPath,fileList(iFile).name],'SMF');
        CurrentEnsemble(iR,:,:) = I.Current;
        SMFEnsemble(iR,:,:) = SMF.SMF;
        
        iR = iR+1;
    end
    
end
load([currentPath,fileList(iFile).name]);


%%


%% CDF

colors='rgbckym'
markers = 'ds*o'
iE = 1;
Estrength_vec =  logspace(-5,3,1000);

sigma_vec = linspace(0,10,50);
iS =1;30;


figure(1);
axes('FontSize',24);
hold on;
grid on;
c=1;
for iE = 1;[201,601,801]
    
    I = log10(abs(CurrentEnsemble(:,iE,iS)));
    %     I = (I- mean(I))/std(I);
    x = linspace(min(I),max(I),length(I));
    
    plot(sort(I),log10((1:length(I))/length(I)),[markers(c),colors(c)],'MarkerSize',8);
    
    c=c+1;
    
end

c=1;
for iE =1; [201,601,801]
    
    I = log10(abs(CurrentEnsemble(:,iE,iS)));
    %     I = (I- mean(I))/std(I);
    x = linspace(min(I),max(I),length(I));
    
    %     plot(x,log10(normcdf(x,mean(I),std(I))),['-',colors(c)],'LineWidth',4);
    plot(x,log10(normcdf(x,mean(I),std(I))),'--k','LineWidth',2);
    plot(x, 0.5*((x-mean(I))/std(I))-log10(2),'LineWidth',2);
    
    c=c+1;
    
end

% title(['\sigma = ',num2str(sigma_vec(iS))]);

legend(['Linear regime          ';...
    'Sinai regime           ';...
    'Saturation             ';...
    'Log-normal distribution';...
    '~x^{1/2}               '],'Location','SouthEast');
xlabel('log10(|I|)');
ylabel('log10(CDF)');
% print(gcf, '-depsc2', 'CDF2');

figure(2);
axes('FontSize',24);
hold on;
grid on;
c=1;
for iE = 1;[201,601,801]
    
    I = log10(abs(CurrentEnsemble(:,iE,iS)));
    %     I = (I- mean(I))/std(I);
    x = linspace(min(I),max(I),length(I));
    plot(sort(I),log10((length(I)-(1:length(I)))/length(I)),[markers(c),colors(c)],'MarkerSize',8);
    
    c=c+1;
    
    
end

% plot(x,log10(1-normcdf(x,0,1)),'-k','LineWidth',4);

% title(['\sigma = ',num2str(sigma_vec(iS))]);

c=1;
for iE = 1;[201,601,801]
    
    I = log10(abs(CurrentEnsemble(:,iE,iS)));
    %     I = (I- mean(I))/std(I);
    x = linspace(min(I),max(I),length(I));
    %     plot(x,log10(1-normcdf(x,mean(I),std(I))),['-',colors(c)],'LineWidth',4);
    plot(x,log10(1-normcdf(x,mean(I),std(I))),'--k','LineWidth',2);
    c=c+1;
    
    
end
xlabel('log10(|I|)');
ylabel('log10(CDF^{-1})');
legend( ['Linear regime          ';...
    'Sinai regime           ';...
    'Saturation             ';...
    'Log-normal distribution'],'Location','SouthWest');
% print(gcf, '-depsc2', 'invCDF2');
%%
figure;
axes('FontSize',24);
hold on;
grid on;
I = CurrentEnsemble(:,iE,iS);

x = linspace(min(I),max(I),length(I));
plot(sort(I),(1:length(I))/length(I),'r','LineWidth',2);
plot(x,normcdf(x,mean(I),std(I)),'LineWidth',2);
xlabel('I');
ylabel('CDF');
legend( ['Data  ';
    'Normal'],'Location','SouthWest');
print(gcf, '-depsc2', 'ISinai_1');

figure;
axes('FontSize',24);
hold on;
grid on;
I = CurrentEnsemble(:,iE,iS);
x = linspace(min(I),max(I),length(I));
plot(sort(I),log10((1:length(I))/length(I)),'r','LineWidth',2);
plot(x,log10(normcdf(x,mean(I),std(I))),'LineWidth',2);
xlabel('I');
ylabel('log10(CDF)');
legend( ['Data  ';
    'Normal'],'Location','SouthWest');
print(gcf, '-depsc2', 'ISinai_2');



figure;
axes('FontSize',24);
hold on;
grid on;
I = CurrentEnsemble(:,iE,iS);
x = linspace(min(I),max(I),length(I));
plot(sort(I),log10(1-(1:length(I))/length(I)),'r','LineWidth',2);
plot(x,log10(1-normcdf(x,mean(I),std(I))),'LineWidth',2);
xlabel('I');
ylabel('log10(CDF^{-1})');
legend( ['Data  ';
    'Normal'],'Location','SouthWest');
print(gcf, '-depsc2', 'ISinai_3');
%%
figure;
axes('FontSize',24);
hold on;
grid on;
I = log10(abs(CurrentEnsemble(:,iE,iS)));
%     I = (I- mean(I))/std(I);
x = linspace(min(I),max(I),length(I));
plot(sort(I),(((1:length(I)))/length(I)),'r','LineWidth',2)
plot(x,(normcdf(x, mean(I),std(I))),'LineWidth',2)
% plot(x, 10.^(0.5*((x-mean(I))/std(I))-log10(2)),'--','LineWidth',4);
plot(x,(0.5*(1-erf((1./x  - mean(1./I))/(sqrt(2)*std(1./I))))),'g--','LineWidth',4)

xlabel('log10|I|');
ylabel('(CDF)');
legend( ['Data            ';
         'Normal          ';...
%          '~x^{1/2}        ';...
         'Inverse Normal  ']);
axis([x(1) x(end) 0 2])
print(gcf, '-depsc2', 'logISinai_1');


figure;
axes('FontSize',24);
hold on;
grid on;
I = log10(abs(CurrentEnsemble(:,iE,iS)));
%     I = (I- mean(I))/std(I);
x = linspace(min(I),max(I),length(I));
plot(sort(I),log10(((1:length(I)))/length(I)),'r','LineWidth',2);
plot(x,log10(normcdf(x, mean(I),std(I))),'LineWidth',2);
% plot(x, (0.5*((x-mean(I))/std(I))-log10(2)),'--','LineWidth',4);
plot(x,(log10(0.5*(1-erf((1./x  - mean(1./I))/(sqrt(2)*std(1./I)))))),'g--','LineWidth',4)



xlabel('log10|I|');
ylabel('log10(CDF)');
legend( ['Data          ';
         'Normal        ';...
%          '~x^{1/2}      ';...
         'Inverse Normal'],'Location','SouthWest');
print(gcf, '-depsc2', 'logISinai_2');
figure;
axes('FontSize',24);
hold on;
grid on;
I = log10(abs(CurrentEnsemble(:,iE,iS)));
%     I = (I- mean(I))/std(I);
x = linspace(min(I),max(I),length(I));
plot(sort(I),log10(1-((1:length(I)))/length(I)),'r','LineWidth',2);
plot(x,log10(1-normcdf(x, mean(I),std(I))),'LineWidth',2);
plot(x,log10(1-0.5*(1-erf((1./x  - mean(1./I))/(sqrt(2)*std(1./I))))),'g--','LineWidth',4)


xlabel('log10|I|');
ylabel('log10(CDF^{-1})');
legend( ['Data          ';
         'Normal        ';
         'Inverse Normal'],'Location','SouthWest');
print(gcf, '-depsc2', 'logISinai_3');
%%
figure;
hold on
SMF = SMFEnsemble(:,201,30);
SMF =(SMF-mean(SMF))/std(SMF)
plot(sort(SMF),log10((1:length(I))/length(I)),linspace(min(SMF),max(SMF),100),log10(normcdf(linspace(min(SMF),max(SMF),100),mean(SMF),std(SMF))));
SMF = SMFEnsemble(:,401,30);
SMF =(SMF-mean(SMF))/std(SMF)
plot(sort(SMF),log10((1:length(I))/length(I)),linspace(min(SMF),max(SMF),100),log10(normcdf(linspace(min(SMF),max(SMF),100),mean(SMF),std(SMF))));
SMF = SMFEnsemble(:,801,30);
SMF =(SMF-mean(SMF))/std(SMF)
plot(sort(SMF),log10((1:length(I))/length(I)),linspace(min(SMF),max(SMF),100),log10(normcdf(linspace(min(SMF),max(SMF),100),mean(SMF),std(SMF))));

