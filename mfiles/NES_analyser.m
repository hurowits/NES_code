%%
figure;
%             hold on
axes('FontSize',24);

subplot(3,2,1:2)
grid on
hold on

plot(SMF_vec/M, (v_drda./D_drda)/2,'b','LineWidth',4);
plot(SMF_vec/M, (V_a./D_a)/2,'-.r','LineWidth',4);
plot(SMF_vec/M, -(velocity./Diffusion)/M/2,'g','LineWidth',4);
plot(SMF_vec/M,1/a_inf*tanh(a_inf*SMF_vec/M/2),'k','LineWidth',4)
plot(SMF_vec/M,1/M*tanh(M*SMF_vec/M/2),'k','LineWidth',4)

set(gca,'FontSize',14)
%     plot(SMF_vec/M, (v_drda./D_drda)/2,'b','LineWidth',4);
%     plot(SMF_vec/M, (V_a./D_a)/2,'-.r','LineWidth',4);
%     plot(SMF_vec/M, -(velocity./Diffusion)/M/2,'g','LineWidth',4);
%     plot(SMF_vec/M,1/a_inf*tanh(a_inf*SMF_vec/M/2),'k','LineWidth',4)
%
%     xlabel('s');
ylabel('v/2D');
legend(['Sample specific     ';...
    'N \rightarrow \infty';...
    'numerics            ';...
    'a_{\infty}          '],'Location','SouthEast');
axis([SMF_vec(1)/M SMF_vec(end)/M 0 1.1/a_inf]);
title(['\sigma=',num2str(sigma),', N=',num2str(M)]);
hold off

subplot(3,2,3)
set(gca,'FontSize',14)

plot(SMF_vec/M,-velocity)
xlabel('s');
ylabel('v');
axtype(2)
subplot(3,2,4)
set(gca,'FontSize',14)

plot(SMF_vec/M,Diffusion)
xlabel('s');
ylabel('D');
axtype(2)



subplot(3,2,5)
set(gca,'FontSize',14)
scatter(2*M*Diffusion,-velocity,20*ones(1,length(SMF_vec)),SMF_vec/M);
ylabel('v');
xlabel('2M*D')

button=1;
s=0;
while (button==1)
    
    
    s_p = sigma*s0(randInd1)+s*ones(1,M);
    w_p = G.*exp(s_p/2);
    w_m = G.*exp(-s_p/2);
    
    subplot(3,2,6)
    set(gca,'FontSize',14)
    
    plot([0,cumsum(log(w_p./w_m))])
    axis([0 M  -M M*sigma])
    xlabel('x');
    ylabel('V(x)');
    title(['s=',num2str(s)]);
    [s,y,button]=ginput(1);
end

