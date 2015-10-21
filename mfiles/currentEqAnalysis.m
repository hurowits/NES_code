sigma = 3;
N = 5;
Delta = 1;
G = genRand(sigma,randperm(N),N+1,'sigma');
energy_vec = linspace(-Delta/2,Delta/2,N);

[w_pl,w_mi] = generateRatesRing (G, energy_vec,N);

num1 = prod(w_pl);
num2 = prod(w_mi);
den=0;
for target = 1:N
    for bond = 1:N
        bond1   = 1 + mod(bond+1,N);
        target1 = 1 + mod(target-1,N);
      den = den + prod(w_pl(bond1:target)) *prod(w_mi(bond:-1:target1));
        [bond1,target,bond,target1]
    end
end