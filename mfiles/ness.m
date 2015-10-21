%[P,P_a,I,er] = ness (w_m, w_p, I_0)
%   computes the non equilibrium steady state distribution and current 
%   in a ring model.
%   Input:
%      w_m, w_p are vectors containing the transition rates. 
%       If I_0 is provided, the calculation is done numerically,
%      otherwise it is done analytically.
%   Output:
%      P - numerical solution
%      P_a - the analytical solution 
%      I - the current
%      er error
function [P,P_a,I,er] = ness (w_m, w_p, I_0)


N = length(w_m);
w_N0 = zeros(1,N);
w_0N = zeros(1,N);

switch nargin
    case 3 %numerical calculation
        
        I = fminsearch(@(I) f_error(I,w_p,w_m),I_0,optimset('TolFun',1e-11,'TolX',1e-11));
        
        P(1) = 1;
        for q = 1:N
            P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
        end
        %
        er=abs(P(end)-P(1)); %error in P
        
        P = P(1:end-1)./sum(P(1:end-1)); %Normalize probabilities
        
        
        I =(P(1)*w_p(1)-P(2)*w_m(1)); %adjust current according to normalization
        P_a=P;
%         er=0;
    case 2 %analytical calculation
        for l = 1:N
            
            V1 = [0 cumsum(log(w_m./w_p))];                 %V1(n)=V(n->0), V1(end) = -SMF
            V2 = [fliplr(cumsum(log(fliplr(w_p./w_m)))) 0]; %V2(n)=V(n->N), V2(1) = SMF
            
            w_N0(l) = 1./sum(exp(V1(1:end-1))./w_p);        %exp(V(n-1->0))/w_{n,n-1}
            w_0N(l) = 1./sum(exp(V2(2:end))./w_m);          %exp(V(n->N))/w_{n-1,n}
            
            P_a(l)=sum(1./w_p.*exp(V1(2:end)));

            
            w_p = circshift(w_p,[0,1]);
            w_m = circshift(w_m,[0,1]);
            
            
            
        end
        P_a=P_a/sum(P_a);
        
        I = 1./sum(1./(w_N0-w_0N));
        P = I./(w_N0-w_0N);
        
        er=sqrt(sum(abs(P-P_a).^2)); %error in P

end

end