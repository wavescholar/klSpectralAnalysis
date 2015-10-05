function [Affty,Lapl,Marov,Kernel] = buildLatentLoop(A,logpow,sizeIm,noScales)
  

  [sA,D,sL,sM] = normalizeAffty(A);
  u0 = D/sum(D);
  
  Affty  = cell(noScales,1);
  Lapl   = cell(noScales,1);
  Markov = cell(noScales,1);
  if noScales > 1
    Kernel = cell(noScales-1,1);
  else
    Kernel = [];
  end
  
  
  logpow = 2;
  [L1,A1,R1,st1,W1,rbinNhbr1,K1,sId1,sMp1] = buildLatent(A,logpow,sizeIm);
  
