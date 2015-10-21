function f = myfun(a_s,s,vd_ratio)
f = 2/a_s * tanh(a_s*s/2) - vd_ratio;
end