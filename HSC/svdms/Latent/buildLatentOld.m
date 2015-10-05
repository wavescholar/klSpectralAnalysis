function [sL,Ar,R,st,W,rbinNhbr,K,selectId,sMp] = buildLatentOld(A,logpow,sizeIm) 
%
%
% 
  
  %=== normalize affinities and generate side info ===
  [sA,D,sL,sM] = normalizeAffty(A);
  sqrtD    = sparse(diag(D.^0.5));
  sqrtDinv = sparse(diag(D .^ -0.5));
  
  %=== kernel selection ===
  [K,selectId,sMp] = latentKernels(sM,sqrtD,sqrtDinv,logpow,sizeIm);
  
  %=== set up latent space probabilities (stationary + markov)
  u0 = D/sum(D) ;
  % old code
  %[Ar,R,st,W,rbinNhbr] = latentProbs(K,D/sum(D));
  % latent Probs now returns K
  [Ar,R,st,W,rbinNhbr,K] = latentProbs(K,D/sum(D));  
  
