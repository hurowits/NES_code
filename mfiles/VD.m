%[v,D] = VD(w_p,w_m)
%   VD Calculates the velocity and diffusion coefficient
%   of a peridoic one dimensional hopping model with transition rates
%   w_{n+1,n}= w_p and w_{n,n+1}= w_m
function [v,D] = VD(w_p,w_m)
M = length(w_p);
W = -diag(w_p + w_m([M,1:M-1]),0) + ...
    diag(w_m(1:M-1),1)+diag(w_p(1:M-1),-1);
W(1,M) = w_p(M);
W(M,1) = w_m(M);



delta_k = .01;

if (M==2)
    
    W(1,2) = w_p(2) + w_m(1);
    W(2,1) = w_m(2) + w_p(1);
    %                 [g0,lambda0,C0] = charFn(W,0);
    lambda0 = charFn(W,0);
    
    W(1,2) = w_p(2)*exp(1i*delta_k*M)  + w_m(1);
    W(2,1) = w_m(2)*exp(-1i*delta_k*M) + w_p(1);
    lambda1 = charFn(W,0);
    %               [g1,lambda1,C1] = charFn(W,0);
    
    W(1,2) = w_p(2)*exp(1i*2*delta_k*M)  + w_m(1);
    W(2,1) = w_m(2)*exp(-1i*2*delta_k*M) + w_p(1);
    %               [g2,lambda2,C2] = charFn(W,0);
    lambda2 = charFn(W,0);
else
    
    lambda1_m = charFn(W,-delta_k/2);
    lambda1_p = charFn(W,delta_k/2);
    lambda2_m = charFn(W,-delta_k);
    lambda2_p = charFn(W,delta_k);
    lambda0 = 0;
end
% lambda0=0;
%dCdk(q) = (C1-C0)/delta_k;
dlambdadk = (lambda1_p - lambda1_m)/delta_k;
%
v = 1i*dlambdadk;
D = -0.5*(lambda2_p + lambda2_m )/delta_k^2 ;


%
%             v_drda(q,iSig) = (1- mean(w_m./w_p))./(mean(1./w_p));
%             D_drda(q,iSig) = 0.5*(mean(1./w_p)).^(-3)*...
%                 (1-(mean(w_m./w_p)).^2)./...
%                 (1-mean((w_m./w_p).^2))*...
%                 (mean(1./w_p.^2)*(1-mean(w_m./w_p))+...
%                 2*mean(w_m./w_p.^2)*mean(1./w_p));

end
