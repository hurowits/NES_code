%f_error(I,w_p,w_m)
%   the function x = P(N+1)-P(1) is to be minimized by fminsearch
%   with respect to I.
%   w_p, w_m are forward and backward transition rates as defined in NER,NEF
function x = f_error (I,w_p,w_m)
P = 1;
for q=1:length(w_p)
%     P(q+1) = (w_p(q)*P(q) - I) / w_m(q);
    P = (w_p(q)*P - I) / w_m(q);
end

x = abs(1-P);

end