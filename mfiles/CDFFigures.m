function [mu,sigma] =  CDFFigures(data,fileName,flag,flagPrint)

if flag==1
    x=linspace(min(data),max(data),100);
    
    figure;
    axes('FontSize',24);
    hold on;
    grid on;
    %     plot(x,(normcdf(x,mean(data),std(data))),'b','LineWidth',4);
    plot(sort(data),((1:length(data))/length(data)),'--k','LineWidth',4);
    
    legend([%'Normal';...
        'Data  '],'Location','SouthEast');
    xlabel(fileName);
    ylabel('CDF');
    axis([x(1) x(end) 0 1])
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'1']);
    end
    %-----------------------------
    figure;
    axes('FontSize',24);
    hold on;
    grid on;
    %     plot(x,log10(normcdf(x,mean(data),std(data))),'b','LineWidth',4);
    plot(sort(data),log10((1:length(data))/length(data)),'--k','LineWidth',4);
    
    legend([%'Normal        ';...
        'Data          '],'Location','SouthEast');
    
    xlabel(fileName);
    ylabel('log10(CDF)');
    axis([x(1) x(end) -3 0])
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'2']);
    end
    %-----------------------------
    %     p = polyfit(sort(data)',log10((length(data):-1:1)/length(data)),2);
    
    figure;
    axes('FontSize',24);
    hold on;
    grid on
    %     plot(x,polyval(p,x),'-b','LineWidth',4);
    %     plot(x,log10(1-normcdf(x,mean(data),std(data))),'-b','LineWidth',4);
    plot(sort(data),log10((length(data):-1:1)/length(data)),'--k','LineWidth',4);
    
    legend([%'Normal';...
        'Data  ']);
    xlabel(fileName);
    ylabel('log10(CDF^{-1})');
    axis([x(1) x(end) -3 0]);
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'3']);
    end
else
    x=linspace(min(data),max(data),100);
    x=linspace(0,max(data),100);
    
    %-----------------------------
    figure;
    axes('FontSize',24);
    hold on;grid on
    p = polyfit(sort(data)',log10((length(data):-1:1)/length(data)),2);
%         p = polyfit(sort(data)',log10((1:length(data))/length(data)),2);
    
    sigma = sqrt(-1/2/(p(1)*log(10)));
    mu = -p(2)/2/(p(1))*log(10);
    
    plot(x,polyval(p,x),'b','LineWidth',4);
    plot(x,log10(1-normcdf(x,mu,sigma)) ,'r','LineWidth',4);

%     plot(x,log10(1-normcdf(x,mu,sigma)) +0.5* log10(sigma^2/2/pi) - mu^2/2/sigma^2/log(10),'--b','LineWidth',4);

    plot(sort(data),log10((length(data):-1:1)/length(data)),'--k','LineWidth',4);
    
    legend(['Normal O(x^2)   ';...
            'Normal CDF (erf)';...
%             '                ';...
            'Data            '],'Location','SouthWest');
    xlabel(fileName);
    ylabel('log10(CDF^{-1})');
    axis([x(1) x(end) -3 0.5]);
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'3']);
    end
    %------------------------------
    figure;
    axes('FontSize',24);
    hold on;grid on
%     plot(x,10.^polyval(p2,x),'b','LineWidth',4);
    plot(x,1-10.^(polyval(p,x)),'b','LineWidth',4);
%         plot(x,(normcdf(x,mean(data),std(data))),'b','LineWidth',4);
    plot(x,normcdf(x,mu,sigma),'r','LineWidth',4);
    plot(sort(data),((1:length(data))/length(data)),'--k','LineWidth',4);
    legend([%'Normal (small I)~O(x^3)';...
        'Normal O(x^2)   ';...
        'Normal CDF (erf)';...
        'Data            '],'Location','SouthEast');
    xlabel(fileName);
    ylabel('CDF');
    axis([x(1) x(end) 0 1])
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'1']);
    end
     %-----------------------------
    figure;
    axes('FontSize',24);
    hold on;grid on
    p2 = polyfit(sort(data(data<-2))',log10((1:length(data(data<-2)))/length(data(data<-2))),3);
%     p2 = polyfit(sort(data)',log10((1:length(data))/length(data)),3);

    
    plot(x,log10(1-10.^polyval(p,x)),'b','LineWidth',4);
    plot(x,log10(normcdf(x,mu,sigma)),'r','LineWidth',4);
    plot(sort(data),log10((1:length(data))/length(data)),'--k','LineWidth',4);
    legend(['Normal O(x^2)   ';...
            'Normal CDF (erf)';...
            'Data            '],'Location','SouthEast');
    
    xlabel(fileName);
    ylabel('log10(CDF)');
    axis([x(1) x(end) -3 0])
    if(flagPrint==1)
        print(gcf, '-depsc2', [fileName,'2']);
    end
end