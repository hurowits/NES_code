%lambda_0 = CHARFN 
%   Calculates the characteristic function of the distribution
%   of a random walk described by the rate equation 
%   dp/dt = Wp
%   Input:
%       W is the transition rate matrix, must be a square matrix.
%       k is the wave number at which the characteristic function is evaluated.
%
%   Output:
%       lambda0 - eigenvalue corresponding to steady state.

% function [ g,lambda_0,C_0 ] = charFn( W, k, t)
% function [ g,lambda_0,C_0 ] = charFn( W, k)
function [lambda_0 ] = charFn( W, k)

M = length(W);
W(1,M) =W(1,M)*exp(-1i*k);
W(M,1) = W(M,1)*exp(1i*k);

lambda = eig(W); %Cols of R are right eigenvectors
% lambda = diag(lambda);
% L = inv(R); %rows of L are left eigenvectors

% M=length(R);
% psi_0= [1;zeros(M-1,1)]; %Initial wave packet

% g=0;


% 
% for nu=M:M
%     g = g + ones(1,M)*R(:,nu) * L(nu,:)*psi_0*exp(lambda(nu)*t);
% end

%Largest eigenvalue and corresponding C
% lambda_0=lambda(M);
[lambda_0,I]=min(abs(lambda));
lambda_0=lambda(I);
% C_0 = ones(1,M)*R(:,I) * L(I,:)*psi_0;



end

