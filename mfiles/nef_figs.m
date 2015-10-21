% --------------
% SMF vs. tau
% --------------
[Gsrt,ind]=sort(G,'Descend');
SMF_d=-cumsum(Delta_n(ind))*beta; %SMF vs. nu for sigma=N

figure;
axes('FontSize',24);
hold on
plot(linspace(0,1,1000),SMF_d,'k','LineWidth',3)
plot(tau_vec,squeeze(SMFEnsemble(1,:,100)),'r','LineWidth',4)
plot(tau_vec,squeeze(SMFEnsemble(1,:,20)),'b','LineWidth',8)
plot(tau_vec,squeeze(SMFEnsemble(1,:,5)),'-.g','LineWidth',8)

% legend(['Random Walk   ';...
%     'SMF, \sigma=50';...
%     'SMF, \sigma=10']);%,'Location','SouthEast');
grid on
% axis([log10(Estrength_vec(1)^2) log10(Estrength_vec(end)^2) -2.5 6])
axis([-0.14 1.2 min(SMF_d(:)) max(SMF_d(:))])
xlabel('ln\nu  [scaled]');
ylabel('SMF')
% print(gcf, '-depsc2', 'SMF_RW.eps');

%%
% -------------------------------------------------------------------------
% I(tau;sigma) image
% -------------------------------------------------------------------------
figure;
axes('FontSize',24);
I=squeeze(CurrentEnsemble(1,:,:));
I(1)=0.02;I(2)=-0.02;
imagesc(tau_vec,sigma_vec(1:90),sign(I(:,1:90)').*(abs(I(:,1:90)))'.^(1/4));
xlabel('ln\nu  [scaled]');
ylabel('\sigma');
% print(gcf, '-depsc2', 'I_sig_tau.eps');

%%
% -------------------------------------------------------------------------
% I(tau;sigma) plots
% -------------------------------------------------------------------------

figure;
axes('FontSize',24);
hold on
plot(log10(nu(2,:).^2*max(G))/sigma_vec(2),-squeeze(CurrentEnsemble(1,:,2))/4,'b','LineWidth',4);
plot(log10(nu(2,:).^2*max(G))/sigma_vec(2),-squeeze(CurrentEnsemble(1,:,1)),'g','LineWidth',8);
grid on
% axis([log10(Estrength_vec(1)^2) log10(Estrength_vec(end)^2) -min(CurrentEnsemble)*1.1 -1.1*max(CurrentEnsemble)])
% axis([-0.2 1.2 -2.5e-2 2e-2])
xlabel('log\nu [scaled]');
ylabel('I');
axis([-0.1 1.1 -5.6e-3 6.4e-3])
legend(['\sigma=50';'\sigma=10'],'Location','SouthWest')
% axtype(1)
% print(gcf, '-depsc2', 'IvsNu.eps');

%%
% -------------------------------------------------------------------------
% SMF and Current log-log
% -------------------------------------------------------------------------
iSigma=1;
SMF = SMFEnsemble(1,:);
Current=CurrentEnsemble(1,:);
figure;
%hold on
axes('FontSize',24);
loglog(Estrength_vec.^2,abs(SMF),'-r', 'LineWidth',4);hold on
loglog(Estrength_vec.^2,abs(Current),'--b', 'LineWidth',4);
% loglog(Estrength_vec.^2,abs(CurrentAnalytical(1,:,2)),':r', 'LineWidth',2);

%     loglog(Estrength_vec.^2,abs(Current_Estimate),'b--', 'LineWidth',2);hold on
% loglog(Estrength_vec.^2,abs(1./w_harm(:,iSigma)),'--r', 'LineWidth',2);
% loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_0(:,iSigma)),'--r', 'LineWidth',2);
% loglog(repmat(Estrength_vec'.^2,1,1),abs(SMF_inf(:,iSigma)),'--r', 'LineWidth',2);
SinaiLimit = sqrt(N)*Delta/2;
% I_typ = std(G)/sqrt(N) *Delta*beta*Estrength_vec.^2;
I_0 = -beta/N*dot(G(iSigma,:),Delta_n)*Estrength_vec.^2;
I_inf = beta/N*dot(1./G(iSigma,:),Delta_n)*ones(length(Estrength_vec),1)/mean(1./G(iSigma,:));
loglog(Estrength_vec.^2,abs(I_0),':m', 'LineWidth',8);
loglog(Estrength_vec.^2,abs(I_inf),':m','LineWidth',8);
% loglog(Estrength_vec.^2, sqrt(N)*Delta/2.*ones(length(Estrength_vec),1),'-k', 'LineWidth',4);

% loglog(Estrength_vec.^2,2/beta*ones(size(Estrength_vec)),'--m','LineWidth',4);
% loglog(Estrength_vec.^2,abs(I_typ),'--g');
% axis([Estrength_vec(1).^2 Estrength_vec(end).^2
% min(abs(Current(:,iSigma))), 10*SinaiLimit]);
% axis([1e-5 1e5 1e-5 1e1]);
currentLims = linspace(min(abs(Current(1,:))),max(abs(SMF(1,:))));
axis([1e-10 Estrength_vec(end).^2 currentLims(1), 10*SinaiLimit]);

loglog(ones(length(currentLims),1).*1./min(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',4);
loglog(ones(length(currentLims),1).*1./max(G(iSigma,:),[],2),currentLims,'k--', 'LineWidth',4);
%     loglog(ones(length(currentLims),1).*1./mean(G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);
%     loglog(ones(length(currentLims),1).*mean(1./G(iSigma,:),2),currentLims,'k:', 'LineWidth',2);

legend(['|SMF|    ';...
    '|Current|';...
    'Eq. 64   ']);
xlabel('\nu');
ylabel('|Current| and |SMF|');
grid on
hold off
% print(gcf, '-depsc2', 'IvsNu.eps');

%
%% Scatter Plots
iE=1;
[SMF,index] = sort(SMFEnsemble(:,iE));
I = CurrentEnsemble(index,iE);
% Select cross section of I vs SMF
smf1 = 2;1e-3;
smf2 = 2.1;1.1e-3;
figure;
axes('FontSize',24);
plot(SMF,(-(CurrentEnsemble(index,iE))),'.')
line([smf1 smf1],[-0.015 0.015],'color','k');
line([smf2 smf2],[-0.015 0.015],'color','k');
% line([-2 -2],[0 0.015],'color','k');
% line([-2.3 -2.3],[0 0.015],'color','k');
axis([-10 10 -0.015 0.015])
grid on
xlabel('SMF')
ylabel('I');
% print(gcf, '-depsc2', 'IvsSMF2.eps');

% P(I,Barrier;SMF)
BarrierSrt = Barrier(index);
b=find((abs(SMF)<smf2).*(abs(SMF)>smf1));

figure;
axes('FontSize',24);
plot(sort(abs(I(b))),1:length(b),'LineWidth',4);
grid on
xlabel('I(2<|SMF|<2.1)');
ylabel('CDF');
%%
% -------------------------------------------------------------------------
% Plots of:
% - Histogram of first max of I
% - Histogram of Barriers
% - Analytical line for Barriers
% -------------------------------------------------------------------------

x = linspace(-15,2,100);
dx=x(2)-x(1);
figure;
axes('FontSize',24);
hold on

% Barrier1=Barrier/std(Barrier);
Barrier1=Barrier/std(Barrier);
hB=hist(-Barrier1 ,x);
[mB,mBind]=max(hB);


I1 = (log(abs(I(b))) - mean(log(abs(I(b)))))/std(log(abs(I(b))));
hI1=hist(I1,x);
[mI,mIind]=max(hI1);
% I1 = I1 - (Barrier1(mBind)+I1(mIind));
% h=hist(I1,x);

plot(x- (-x(mBind)+x(mIind)),hI1/sum(hI1)/dx,'or','MarkerSize',12,'LineWidth',2)
% plot(x,h/sum(h)/dx,'-or','MarkerSize',8)
maxI1 = (log(abs(maxI)) - mean(log(abs(maxI))))/std(log(abs(maxI)));
hMaxI=hist((maxI1),x);
[mmI,mmIind]=max(hMaxI);
% maxI1 = maxI1 - (Barrier1(mBind)+maxI1(mmIind));
% h=hist(maxI1,x);

plot(x- (-x(mBind)+x(mmIind)),hMaxI/sum(hMaxI)/dx,'db','MarkerSize',12,'LineWidth',2)
% plot(x,h/sum(h)/dx,'-db','MarkerSize',8)

plot(x,hB/sum(hB)/dx,'*g','MarkerSize',12,'LineWidth',2);

sigmaU=sqrt(N*(std(Delta_n)*beta)^2*log10(max(G)*Estrength_vec(iE)^2)/(sigma));

a=pi^2*sigmaU^2/2;

x2=-x;
f=zeros(1,length(x2));
for n=1:2:101
    f = f+ 2*a./x2.^3.*exp(-a*n^2./x2.^2).*(2*a*n^2./x2.^2-1)*8/pi^2;
    
end

plot(x,f,'--k','LineWidth',4)


legend(['P(I_{max})';...
    'P(I;SMF)  ';...
    'P(B)      ';...
    'Theory    '],'Location','NorthWest');
grid on
xlabel('ln(|I|)  [scaled]');%(2<|SMF|<2.3)');
ylabel('Probability density');
axis([-10 0 0 0.6 ])
% print(gcf, '-depsc2', 'PlnI.eps');


%%
figure;
axes('FontSize',24);
hold on
plot(sort(BarrierSrt(b)),1:length(b),'LineWidth',4)
% plot(sort(-log(abs(I(b)))),1:length(b),'r','LineWidth',4)
% hist(abs(I(b)),100);
grid on
xlabel('B');%(2<|SMF|<2.3)');
ylabel('CDF(B)');

figure;
axes('FontSize',24);
[h,bins]=hist((-(BarrierSrt(b))),100);
plot(bins,h/sum(h),'-r','LineWidth',4)

% hist(abs(I(b)),100);
grid on
xlabel('-B');%(2<|SMF|<2.3)');
ylabel('Prob(B)');
axis([-6 0 0 0.035])
% print(gcf, '-depsc2', 'PBarrier.eps');


figure;
axes('FontSize',24);
hold on
plot(-BarrierSrt(b)/2,log(abs(I(b))),'.')
plot([-2.5 -1], [-7 -5.5], 'k', 'LineWidth',4)
grid on
xlabel('-B');
ylabel('ln(|I|)');%(2<|SMF|<2.3)');
% print(gcf, '-depsc2', 'BvsLnI.eps');




figure;
hist(-BarrierSrt(b),100)
minVsrt=minV(index);
maxVsrt=maxV(index);
[hmin,binsmin]=hist(minVsrt(b),100);
[hmax,binsmax]=hist(maxVsrt(b),100);

figure;
axes('FontSize',24);

hold on
plot(maxVsrt(b),minVsrt(b),'g.');
plot(10*hmin/sum(hmin),binsmin,'b','LineWidth',2)
plot(binsmax,-10*hmax/sum(hmax),'r','LineWidth',2)
grid on
xlabel('u^{+}');
ylabel('u^{-}');
legend(['(u^+,u^-)';...
    'p(u^-)   ';...
    'p(u^+)   '], 'Location','SouthEast')
axis([0 4 -4 0])
% print(gcf, '-depsc2', 'upum.eps');

%%
sigmaU=sqrt(N*(std(Delta_n/2)*beta)^2*log10(max(G)*Estrength_vec(iE)^2)/(sigma));
figure;
axes('FontSize',24);
hold on
x = (pi*sigmaU)^2/2./(sort(BarrierSrt(b))).^2;
plot(-x,log((1:length(BarrierSrt(b)))/length(b))-pi/2,'.','MarkerSize',24);
plot([-10 0],[-10 0],'k','LineWidth',4)
legend(['Data';...
    'y=x '],'Location','NorthWest');
xlabel('-[\pi^2\sigma_U^2/8]/B^2')
ylabel('ln [Prob(Barrier<B)]')
legend(['Data';...
    'y=x '],'Location','Best');
grid on
% print(gcf, '-depsc2', 'CDFB.eps');
%
Bf = sort(BarrierSrt(b));
figure;
axes('FontSize',30);hold on
B = linspace(0.0001,15,50);
h=hist(BarrierSrt(b),B);
plot(B,h/sum(h),'-r','LineWidth',4)
PB= pi^2*sigmaU^2*exp(-pi^2*sigmaU^2/2./B.^2)./B.^3*(B(2)-B(1));
% PB=PB/sum(PB);
plot(B, PB,'.','MarkerSize',20)
plot([Bf(10);Bf(10)],[0;0.1],'--k','LineWidth',2)
plot([Bf(3774);Bf(3774)],[0;0.1],'k','LineWidth',2)
%
% B2 = linspace(0.0001,15,25);
% h2=hist(-log(abs(maxI(1:1306))),B2);
% plot(B2,h2/sum(h2))
grid on

xlabel('B');%(2<|SMF|<2.3)');
ylabel('Prob(B)');
% print(gcf, '-depsc2', 'PB2.eps');

% axis([-6 0 0 0.035])
%%
figure;
hist(abs(I(b)),100);
xlabel('|I|');

figure;
hist(1./abs(I(b)),100);
xlabel('1/|I|');

figure;
axes('FontSize',24);
hist(-log(abs(I(b))),100)
xlabel('-log(|I|)');
ylabel('P(log|I| | SMF)')
grid on
% print(gcf, '-depsc2', 'PlogI');

figure;
axes('FontSize',24);
hist(Barrier(:),100);
xlabel('B');
ylabel('P(B)');
grid on

% print(gcf, '-depsc2', 'PmiB');
% plot(x,log10(numRealizations*(1-0.5*(erf((mu+x)/sqrt(2*sigma^2))-erf((mu-x)/sqrt(2*sigma^2))))),'--r','LineWidth',4);


%%
I=((CurrentEnsemble(index,iE)));
figure;
plot(SMF,I)
p = polyfit(SMF,I,1);
%remove mean
I = I - polyval(p,SMF);
%remove variance
sigma=sqrt(conv(I.^2,ones(1,100),'same')/100);
I = (I)./sigma;
figure;plot(SMF, abs(I),'.');
%%


%% CDF

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort(abs(I)),1:numRealizations,'.','LineWidth',2);
xlabel('|I|');
ylabel('CDF');

figure;
axes('FontSize',24);
hold on;
grid on
CurrentNormalized = (CurrentEnsemble - repmat(mean(CurrentEnsemble),[numRealizations 1]))./(repmat(std(CurrentEnsemble),[numRealizations 1]));
plot(sort(abs(CurrentNormalized)),1:numRealizations,'LineWidth',2);
x= (linspace(min(abs(CurrentNormalized(:))),max(abs(CurrentNormalized(:))),1000));

mu = mean(CurrentNormalized);
sigma = std(CurrentNormalized);
plot(x,(numRealizations*(0.5*(erf((mu(1)+x)/sqrt(2*sigma(1)^2))-erf((mu(1)-x)/sqrt(2*sigma(1)^2))))),':b','LineWidth',4);
plot(x,(numRealizations*(0.5*(erf((mu(2)+x)/sqrt(2*sigma(2)^2))-erf((mu(2)-x)/sqrt(2*sigma(2)^2))))),':g','LineWidth',4);
plot(x,(numRealizations*(0.5*(erf((mu(3)+x)/sqrt(2*sigma(3)^2))-erf((mu(3)-x)/sqrt(2*sigma(3)^2))))),':r','LineWidth',4);
%
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,1)),std(CurrentNormalized(:,1))),'b','LineWidth',2)
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,2)),std(CurrentNormalized(:,2))),'g','LineWidth',2)
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,3)),std(CurrentNormalized(:,3))),'r','LineWidth',2)
xlabel('|I|');
ylabel('CDF');
%% Inverse CDF

figure;
axes('FontSize',24);
hold on;
grid on
plot(sort(abs(CurrentEnsemble)),numRealizations:-1:1,'.','LineWidth',2);
xlabel('|I|');
ylabel('CDF^{-1}');

figure;
axes('FontSize',24);
hold on;
grid on
CurrentNormalized = (CurrentEnsemble - repmat(mean(CurrentEnsemble),[numRealizations 1]))./(repmat(std(CurrentEnsemble),[numRealizations 1]));
plot(sort(abs(CurrentNormalized)),numRealizations:-1:1,'LineWidth',2);
x= (linspace(min(abs(CurrentNormalized(:))),max(abs(CurrentNormalized(:))),1000));

mu = mean(CurrentNormalized);
sigma = std(CurrentNormalized);
plot(x,(numRealizations*(1-0.5*(erf((mu(1)+x)/sqrt(2*sigma(1)^2))-erf((mu(1)-x)/sqrt(2*sigma(1)^2))))),':b','LineWidth',4);
plot(x,(numRealizations*(1-0.5*(erf((mu(2)+x)/sqrt(2*sigma(2)^2))-erf((mu(2)-x)/sqrt(2*sigma(2)^2))))),':g','LineWidth',4);
plot(x,(numRealizations*(1-0.5*(erf((mu(3)+x)/sqrt(2*sigma(3)^2))-erf((mu(3)-x)/sqrt(2*sigma(3)^2))))),':r','LineWidth',4);
%
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,1)),std(CurrentNormalized(:,1))),'b','LineWidth',2)
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,2)),std(CurrentNormalized(:,2))),'g','LineWidth',2)
%   plot(x,numRealizations*normcdf(x,mean(CurrentNormalized(:,3)),std(CurrentNormalized(:,3))),'r','LineWidth',2)

xlabel('|I|');
ylabel('CDF^{-1}');

% plot(sort(abs(CurrentEnsemble2)),log10(numRealizations:-1:1),'.k');

%             figure;
%             axes('FontSize',24);
%             hold on;
%             grid on
%             plot(1:N, V1(1:end-1),'b','LineWidth',2);
%             plot(1:N,squeeze(P_analytical(1,1,end:-1:1))*N,'g','LineWidth',2);
%             plot(1:N,N*exp(-V1(1:end-1))./sum(exp(-V1(1:end-1))),'r','LineWidth',2);
%
%             legend(['V(x)       ';...
%                     ' N*P(x)    ';...
%                     '~exp(-V(x))'],'location','NorthWest')
%             xlabel('x');
% %             print(gcf, '-depsc2', 'PvsV2');
%%
% Vflat = V1-V1(end)*(0:N)/(N);
% Vflat = Vflat(2:end) -0.5*log(w_p./w_m);
% Vflat =  Vflat(1:end-1);
% indshift =
iE=1;295;
figure;
axes('FontSize',24);
hold on;
grid on
% plot(1:N, V1(2:end),'k','LineWidth',2);
plot(1:N,fftshift(squeeze(P_Numerical(1,iE,:))),'k','LineWidth',4);

plot(1:N,fftshift(exp(-U(iE,:))./sum(exp(-U(iE,:)))),'--r','LineWidth',4);


w = sqrt(w_m.*w_p);
f=zeros(1,N);
f2=zeros(1,N);
w_g=zeros(1,N);

z=zeros(1,N);
xx=0:N-1;
SMF=SMFEnsemble(iE);
for q=1:N
    f(q) = sum(exp(U(iE,:) -SMF/N*xx)./w);
    %     f2(q) = sum(exp(Vflat +sign(SMF)*xx)./w);
    z(q) = sum(exp(U(iE,:) -SMF/N*xx));
    w_g(q) = f(q)/z(q);
    w=circshift(w,[0 -1]);
    xx=circshift(xx,[0 -1]);
    %     Vflat = circshift(Vflat,[0 -1]);
end
plot(1:N, fftshift(U(iE,:)/N - 4e-3),'--c','LineWidth',8);
plot(1:N,fftshift(fftshift(-log(z/sum(z)))/N)-10e-3,'b','LineWidth',2);
% plot(1:N,fftshift(fftshift(-log(f/sum(f)))/N)-12e-3,':r','LineWidth',2);
% plot(1:N,fftshift(U(iE,:)+log(z))/N-12e-3)
%   plot(1:N,(1./w)/sum(1./w)*N,'y','LineWidth',2);

% plot(1:N,log(z),'y','LineWidth',2);
%  plot(1:N,w_g/sum(w_g)*N,'--k','LineWidth',2);
% plot(1:N,exp(-Vflat).*f/sum(exp(-Vflat).*f)*N,':m','LineWidth',2);
% plot(1:N,1./(N*exp(-Vflat(1:end))./sum(exp(-Vflat(1:end))) ./  (squeeze(P_analytical(1,1,end:-1:1))*N)'),':k','LineWidth',2)
% plot(1:N,exp(-Vflat).*f2/sum(exp(-Vflat).*f2)*N,':c','LineWidth',2);
% plot(1:N,N*P,'--c','LineWidth',2);
hl1=legend([
    %'V(x)             ';...
    %'V(x) - SMF x/N';...
    'P(x)               ';...
    'e^{-U(x)}          ';...
    %     '(1/w(x))_{\epsilon}';...
    'U(x)               ';...
    'U_{\epsilon}(x)    ']);
%     '(1/w)_{\epsilon}   ']);
%     '~e^{-U(x)}<<1/w>>';...
%     'ref              '],...
%
set(hl1,'Interpreter','TeX')
xlabel('x');
ylabel('Prob(x), U(x)')
% axis([1 N -4.1 4.1]);
% print(gcf, '-depsc2', 'PvsV4');