function [Ar,R,st,W,rbinNhbr,K] = latentProbs(K,u0)

  [Ar,R,st,W,K] = setLatent(K,u0,1) ; % normalizes Ar by the median
  rbinNhbr = Ar > 0;
  rbinNhbr = rbinNhbr - eye(size(rbinNhbr)) ;
