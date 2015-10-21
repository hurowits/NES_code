function [P,I,SMF,er] = ness (w_m, w_p, I_0)

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
        er(iE)=abs(P(end)-P(1)); %error in P
        
        P = P(1:end-1)./sum(P(1:end-1)); %Normalize probabilities
        
        
        I =(P(1)*w_p(1)-P(2)*w_m(1))*N ; %adjust current according to normalization
        
    case 2 %analytical calculation
        for l = 1:N
            
            V1 = [0 cumsum(log(w_m./w_p))];                 %V1(n)=V(n->0), V1(end) = -SMF
            V2 = [fliplr(cumsum(log(fliplr(w_p./w_m)))) 0]; %V2(n)=V(n->N), V2(1) = SMF
            
            w_N0(l) = 1./sum(exp(V1(1:end-1))./w_p);        %exp(V(n-1->0))/w_{n,n-1}
            w_0N(l) = 1./sum(exp(V2(2:end))./w_m);          %exp(V(n->N))/w_{n-1,n}
            
            w_p = circshift(w_p,[0,1]);
            w_m = circshift(w_m,[0,1]);
            
        end
        
        I = 1./sum(1./(w_N0-w_0N));
        P = I./(w_N0-w_0N);
end

SMF = sum(log(w_m)-log(w_p));

end