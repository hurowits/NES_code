% [w_pl,w_mi] = generateRatesRing (G, energy_vec)
% Input: 
% G is a vector of driven transitions
% energy_vec is a vector of on-site energies
function [w_pl,w_mi] = generateRatesRing (G, energy_vec,N)
% [w_pl,w_mi] = generateRatesRing (G, energy_vec)
% Input: 
% G is a vector of driven transitions
% energy_vec is a vector of on-site energies
% Output:
% w_pl, w_mi are total clockwise and counter clockwise transitions
w_beta=1;beta=1;
    delta_vec = [diff(energy_vec), energy_vec(N)-energy_vec(1)];
    w_pl = log(G + 2*w_beta./(1+exp(beta*delta_vec)));
    w_mi = log(G + 2*w_beta./(1+exp(-beta*delta_vec)));
   
end