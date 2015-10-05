function binSupp = nonmaxsupp(dLogHalf,u,binNhbr,dLogHalfMaxCut,binCut0)

FALSE = 0; TRUE = 1;

binSupp = zeros(size(binNhbr));
binSupp(:) = 0;
    
    dLogHalf = reshape(dLogHalf, length(u), length(u)); 
    idx = (dLogHalf(:) < dLogHalfMaxCut) & (~binCut0(:)); 
    if any(idx) 
      tt = find(idx); 
      indIp = ceil(tt/length(u)); 
      indJp = tt - (indIp-1)*length(u); 
      for k=1:length(indIp) 

	suppress = FALSE; 
         
        % Set i0, j0 to be endpoints of edge (indJp(k), indIp(k)) 
        % ordered so that u(i0) > u(j0) 
        if (u(indIp(k)) > u(indJp(k))) 
          i0 = indIp(k); 
          j0 = indJp(k); 
        elseif (u(indJp(k)) > u(indIp(k))) 
          j0 = indIp(k); 
          i0 = indJp(k); 
        else 
          sprintf('Huh? %d %d %f %f\n', indJp(k), indIp(k), u(indIp(k)), dLH) 
          suppress = TRUE; 
        end 
          
        dLH = dLogHalf(j0, i0); 
         
        % Find nhbrs of i0, in current graph with cuts, for which u is 
        % larger than u(i0).  Nonmax suppression occurs if dLogHalf is 
        % still more negative on edges from i0 to these nhbrs. 
        altNhbrs = find(binNhbr(:,i0) & ~binCut0(:, i0) & u(:) > u(i0)); 
        if ~suppress & any(altNhbrs) 
          [mn jmn] = min(dLogHalf(altNhbrs, i0)); jmn = altNhbrs(jmn); 
          if mn < dLH 
            suppress = TRUE; 
          end 
        end 
         
        % Find nhbrs of j0, in current graph with cuts, for which u is 
        % smaller than u(j0).  Nonmax suppression occurs if dLogHalf is 
        % still more negative on edges from j0 to these nhbrs. 
        altNhbrs = find(binNhbr(j0,:) & ~binCut0(j0, :) & (u(j0) > u(:))'); 
        if ~suppress & any(altNhbrs) 
          [mn jmn] = min(dLogHalf(j0, altNhbrs)); jmn = altNhbrs(jmn); 
          if mn < dLH 
            suppress = TRUE; 
          end 
        end 
         
        if suppress 
          binSupp(indJp(k), indIp(k)) = TRUE; 
          binSupp(indIp(k), indJp(k)) = TRUE; 
        end 
      end 
    end  % over idx, all edges with sufficiently small dLogHalf 