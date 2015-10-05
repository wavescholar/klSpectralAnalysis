function [selectId,C] = gramSelect(B,minRes,debug)
% gram selection with random drawing of basis vectors.
%
% see also gramFixed.m

  if (nargin < 3)
    debug = 0;
  end
  
  FALSE = (0 == 1);
  TRUE = ~FALSE;

  %% Params for selecting initial K-means centers 
  %n = size(B);
  n = size(B,2);  
  
  %% D0 center selection: 
  %% Rerun: rand('state', saveRandState); 
  saveRandState = rand('state'); 
  
  %% Normalize cols of sLp;  
  %res = full(sLp * sparse(diag(sum(sLp .^ 2, 1).^-0.5))); 
  res = full(B * sparse(diag(sum(B .^ 2, 1).^-0.5)));   
  
  %% Select initial centers to span space.. 
  j = round(0.5 + n(1) * rand(1)); 
  j = max(j, 1); j = min(j,n(1));
  
  selectId = [j];
  c = res(:, j); 
  c = c/norm(c);
  k = 2;
  done = FALSE; 
  while ~done 
    res = res - c(:,k-1) * (c(:,k-1)' * res);
    nrmRes2 = sum(res .* res, 1); 
    samp = cumsum(nrmRes2 .* (nrmRes2>minRes)); 
    if samp(n(1)) == 0 
      done = TRUE;
      break; 
    end  
    samp = samp/samp(n(1)); 
    r = rand(1); 
    idx = find(samp>=r);    
    
    if any(idx)
      newc = res(:, idx(1));
      newc = newc/norm(newc); 
      c = [c  newc]; 
      selectId = [selectId, idx(1)];
      if (debug) 
	fprintf(2, ' k %d, basis # %d\n', k, idx(1));
      end
      
      k = k+1; 
    else 
      error('Random draw fanned!??!'); 
    end  
  end 
  k = k-1; 
  dm = k;
  C = c;
  %fprintf(2,' Got %d basis elements for power %f\n', k, 2^floor(logpow)); 
  
