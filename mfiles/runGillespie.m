% clear all
sigma =2;
Delta=1;
N=10;
beta=0.1;
w_beta=1;

SMF_vec=[linspace(0,20,15)];%, linspace(5,50,10)];


% SMF_vec=[linspace(0,80,50),80:10:150];%80:10:150;

% SMF_vec=20;
% SMF_vec=0:5;
% for q=51:51+length(SMF_vec)-1
for q=1:length(SMF_vec)
    
    SMF = SMF_vec(q); %SMF per unit cell!
    
    
%     M=1, AAAAA
        M = 1;
        w_p = ones(1,M)*exp(SMF);
        w_m = ones(1,M);
    
        SMF1(q) = sum(log(w_p./w_m))/M;
        [muQ1(q),sigQ1(q),muX1(q),sigX1(q),t_end1(q)] = Gillespie3(w_p,w_m);
    
    % %     %M=2, ABABA
        M = 2;
        a = 1;
        b = 1;
        s1=SMF/M+sigma;
        s2 = SMF/M-sigma;
        A1 = exp(s1/2);
        A2 = exp(-s1/2);
        B1 = exp(s2/2);
        B2 = exp(-s2/2);
%         A1 = a * exp(SMF/3);
%         A2 = 1;
%         B1 = b * exp(SMF*2/3);
%         B2 = 1; 
        w_p = [A1,B1];
        w_m = [A2,B2];
        SMF2(q) = sum(log(w_p./w_m))/M;
        [muQ2(q),sigQ2(q),muX2(q),sigX2(q),t_end2(q)] = Gillespie3(w_p,w_m);
%         s1 = 2*SMF/3;
%         s2 = 2*2*SMF/3;
%         
%           sigma2(q) = 1- 4*a*b * cosh((s1+s2)/4).^2./(a*cosh(s1/2)+b*cosh(s2/2)).^2;
    
%         rat2(q) = 4*tanh(SMF/4)/(1+sigma(q) * tanh(SMF/4)^2);
    
%         mu2(q) = 2*(w_p(1)*w_p(2)-w_m(1)*w_m(2))/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
%         D2(q) = (2*(w_p(1)*w_p(2)+w_m(1)*w_m(2)) - mu2(q)^2)/(w_p(1)+w_p(2)+w_m(1)+w_m(2));
%     
    
    %     %M=10 ABCDEF....
%     M=10;
% %     w_p = [exp(SMF/2),ones(1,M-1)];
% %     w_m = [exp(-SMF/2),ones(1,M-1)];
%     E = rand(M,1);
%     E=E/sum(E)*SMF;
% %     SMF/(M*(M-1))*(0:M-1);
%     w_p=exp(E(randperm(M))/2);
%     w_m=exp(-E(randperm(M))/2);
%     w_p = exp(SMF/(M*(M-1))*(0:M-1));
%     w_m = exp(-SMF/(M*(M-1))*(0:M-1));
% %     
%     [muQ8(q),sigQ8(q),muX8(q),sigX8(q),t_end8(q)] = Gillespie3(w_p,w_m);
    %
    %
    %
    %
    
    SMF
    q
end

%%
SMF2=linspace(0,20,200);

% s1=2*SMF2/3;s2=2*2*SMF2/3; 
s1 = SMF2+sigma;
s2 = SMF2-sigma;
sigma_est = 1- 4*a*b * cosh((s1+s2)/4).^2./(a*cosh(s1/2)+b*cosh(s2/2)).^2;
rat= 4*tanh(SMF2/4)./(1+sigma_est.*tanh(SMF2/4).^2);
figure;
axes('FontSize',24);
hold on;
grid on


plot(SMF_vec,muX1./sigX1.^2,'k*','MarkerSize',10); %AAAA

plot(SMF_vec/M,abs(muX2)./sigX2.^2,'ro','LineWidth',3); %ABAB
plot(SMF2,SMF2/2,'--b','LineWidth',5)


plot(SMF2,tanh(SMF2/2),'k','LineWidth',5);
plot(SMF2/M,tanh(a_inf*SMF2/M/2)/a_inf,'k','LineWidth',4);

plot(SMF2/M,rat/4,'r','LineWidth',2); %ABAB
plot(SMF2/M,tanh(SMF2/2)/2,'k','LineWidth',5);
% plot(abs(SMF_vec)/M/2,abs(mu2)./(D2*2),'r','LineWidth',3); %random chain


% plot(SMF_vec/2, 2*tanh(SMF_vec/4)./(1+tanh(SMF_vec/4).^2))

% plot(abs(SMF_vec),abs(mu3)./(D3*2),'--g','LineWidth',3); %random chain
% plot(abs(SMF_vec),SMF_vec/2,'k','LineWidth',3); %random chain


xlabel('s');
ylabel('v/2D');
legend(['1 site ';...
        '2 sites';...
        'ESR    '],...
    'Location', 'SouthEast');
axis tight;
axis([10^(-1) 10,0, 1.1])
% axtype(1);
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NES/Figs/v_D_sim_new.eps');