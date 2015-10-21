function axtype(option)

% complementary to 'axis'.
% option= 0, 1, 2, 3  for 
% normal,logx,logy,loglog.

% written by Doron Cohen

if option==0;
set(gca,'xscale','linear');
set(gca,'yscale','linear');
end;

if option==1
set(gca,'xscale','log')
set(gca,'yscale','linear')
end

if option==2
set(gca,'xscale','linear')
set(gca,'yscale','log')
end

if option==3
set(gca,'xscale','log')
set(gca,'yscale','log')
end
