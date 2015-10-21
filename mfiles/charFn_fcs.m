% function [ g,lambda_0,C_0 ] = charFn( W, k, t)
% function [ g,lambda_0,C_0 ] = charFn( W, k)
function [P] = charFn_fcs( W, k)
%CHARFN Calculate the characteristic function of the distribution
%of a random walk described by the rate equation
% dp/dt = Wp
% Input:
% W is the transition rate matrix, must be a square matrix.
% k is the wave number at which the characteristic function is evaluated.

M = length(W);
W(1,M) =W(1,M)*exp(-1i*k);
W(M,1) = W(M,1)*exp(1i*k);


if(k==0)
    
    B = [W;ones(1,M)];
    C = zeros(M+1,1);C(M+1)=1;
    [P,R] = linsolve(B,C);
    
else
    B=W;
    C=zeros(M,1);
    [P,R]=linsolve(B,C);
end

end

