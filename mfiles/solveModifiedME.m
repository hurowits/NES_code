
function P = solveModifiedME (w_p, w_m)
%solveModifiedME


N = length(w_m);
lambda_vec=linspace(-pi,pi,N);
time_vec = linspace(0,1e2,100);
for iT = 1:length(time_vec)
    
    for iL = 1:length(lambda_vec)
        lambda = lambda_vec(iL);
        t= time_vec(iT);
        %     W = -diag(w_drv + w_drv([N,1:N-1]) + w_bath_cw + w_bath_ccw([N,1:N-1]),0)+...
        %         diag(w_drv(1:N-1) + w_bath_ccw(1:N-1),1)+...
        %         diag(w_drv(1:N-1) + w_bath_cw(1:N-1),-1);
        %
        %     W = -diag(w_p + w_m([N,1:N-1]),0) + ...
        %         diag(w_m(1:N-1).*  (w_m(1:N-1)./w_p(1:N-1)).^(-lambda),1) + ...
        %         diag(w_p(1:N-1).* (w_p(1:N-1)./w_m(1:N-1)).^(-lambda),-1);
        %
        
        W_lam = -diag(w_p + w_m([N,1:N-1]),0) + ...
            diag(w_m(1:N-1),1) + ...
            diag(w_p(1:N-1),-1);
        W_lam(N,1) = w_m(N);
        W_lam(1,N) = w_p(N);
        W_lam(1,2) = W_lam(1,2)*exp(-1i*lambda);
        W_lam(2,1) = W_lam(2,1)*exp(1i*lambda);
        
        Pbar(iL,iT) = sum(expm(W_lam*t)*[1;zeros(N-1,1)]);
        
%         [V,D] = eig(W_lam);
        
        %
        %
        %
        %     B=[W; ones(1,N)];
        %     C=zeros(N+1,1);
        %     C(N+1)=1;
        %     [psi(iL,:),R]=linsolve(B,C);
        %
        %     p_bar(iL) = sum(psi(iL,:));

    end
    
%     P(:,iT) = ifft(Pbar(:,iT),[],1);
    figure(1);
    plot(lambda_vec,real(fftshift(ifft(Pbar(:,iT)))),'.')
    figure(2);
    plot(lambda_vec,imag(fftshift(ifft(Pbar(:,iT)))),'.')
end
end