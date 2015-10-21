function [mu,sig,Q] = Q_stats (site_t,t1,t2,N)
   
site_t = site_t(t1:t2,:);

for iR = 1:size(site_t,2)
    
    ind1 = find((mod(site_t(:,iR),N)+1)==4);
  
    if(~isempty(ind1) && ind1(end)==size(site_t,1))
        ind1 = ind1(1:end-1);
    end
    if(~isempty(ind1) && ind1(1)==1)
        ind1 = ind1(2:end);
    end
    
    ind2 = find(mod(site_t(ind1-1,iR),N)+1==3);
   
    ind3 = find(mod(site_t(ind1+1,iR),N)+1==3);
    
    Q(iR) = length(ind2)-length(ind3);
    
end
    mu = mean(Q);
    sig = std(Q);
end