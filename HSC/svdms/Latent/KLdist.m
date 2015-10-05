function bts = KLdist(pKnown,qEst)
% probability distribution 
% pKnown : known
% qEst   : estimate
%
  
  small  = 1.0e-08;
  pKnown = pKnown(:);
  qEst   = qEst(:);
  
  if (length(pKnown) ~= length(qEst))
    error('unequal length of the distributions');
  end
  
  if (abs(sum(pKnown)-1.0) > small)
    %fprintf('pKnown does not sum up to 1.0: %f\n',abs(sum(pKnown)-1.0));
    str = sprintf('pKnown does not sum up to 1.0: %f\n',abs(sum(pKnown)-1.0));
    %error(str);
    pKnown = pKnown / sum(pKnown);
  end

  if (abs(sum(qEst)-1.0) > small)
    %fprintf('qEst does not sum up to 1.0: %f\n',abs(sum(qEst)-1.0));
    str = sprintf('qEst does not sum up to 1.0: %f\n',abs(sum(qEst)-1.0));
    %error(str);
    qEst = qEst / sum(qEst);
  end

  
  pp = pKnown .* log2(pKnown+eps);
  pq = pKnown .* log2(qEst+eps);
  bts  = sum(pp - pq);
  
