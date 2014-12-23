function [ns1 ns2] = findfactors(nhh)

if isprime(nhh)
    nhhp = nhh + 1;
else
    nhhp = nhh;
end
fhh = factor(nhhp);
nfhh = floor(length(fhh)/2);
ns1 = 1; ns2 = 1;
for t = 1:nfhh
    if mod(t,2) == 1
        ns1 = fhh(t)*ns1;
        ns2 = fhh(end-t+1)*ns2;
    else
        ns2 = fhh(t)*ns2;
        ns1 = fhh(end-t+1)*ns1;
    end
end
if nfhh ~= length(fhh)/2
    ns1 = ns1*fhh(nfhh+1);
end
a = sort([ns1,ns2]); 
ns1 = a(1);
ns2 = a(2);

% see if just finding sqrt gives better dims:
sq = sqrt(nhh);
if sq - floor(sq) > 0
    sq = floor(sq)+1;
end

if max(ns1,ns2) > sq
    ns1 = sq;
    ns2 = sq;
end
    
    

