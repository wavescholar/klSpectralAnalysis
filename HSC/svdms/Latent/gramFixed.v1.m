function [selectId,C] = gramFixed(B,debug)
% gram selection. but no random drawing.
% the basis vectors are taken in the order
% they are seen in matrix B. 
%
% see also gramSelect.m

  if nargin < 2
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
  %res = full(B * sparse(diag(sum(B .^ 2, 1).^-0.5)));   
  res = zeros(size(B));
  for zz = 1:size(B,2)
    res(:,zz) = B(:,zz)/sqrt(sum(B(:,zz).^2));
  end
  
  
  %% Select initial centers to span space.. 
  initId = 1;
  j = initId;
  
  selectId = [j];
  c = res(:, j); 
  c = c/norm(c);
  k = 2;
  done = FALSE; 
  while ~done 
    res = res - c(:,k-1) * (c(:,k-1)' * res);

    if length(selectId) == size(B,2)
      done = TRUE;
      break; 
    end  

    initId = initId + 1;

    newc = res(:, initId);
    newc = newc/norm(newc); 
    c = [c  newc]; 
    selectId = [selectId, initId];
    if (debug) 
      fprintf(2, ' k %d, basis # %d\n', k, initId);
    end
    
    k = k+1; 
  
  end 
  k = k-1; 
  dm = k;
  C = c;
  %fprintf(2,' Got %d basis elements for power %f\n', k, 2^floor(logpow)); 
  
