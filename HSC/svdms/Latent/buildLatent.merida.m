function [sL,Ar,K,R,st,W,rbinNhbr,selectId,sMd] = buildLatent(A,logpow,sizeIm) 
%
% buildLatent: set up latent space, build kernels, 
%              generate transitional matrix for 
%              the low-dim space.
%
% calls to:
%   normalizeAffty.m 
%   emKernelFitNew.m
%
  %======== STEP 1: SET PARAMETERS ===========

  %=== kernel selection parameters ===
  % find an ordering on the kernels based on max height.
  % increasing this # will increase the total 
  % number of kernels selected
  MAX_HT_FRAC = 0.5;%0.5;
  % length of the intersection 
  INTERSECT_FRAC = 0.4;%0.4;
  MISS_FRAC = 0.75;
  % to threshold small probabilities. look in STEP 4 below.
  SUPP_PROB_FRAC = 0.01;
  
  %% parameter set that works: Apr 23, '03
  %MAX_HT_FRAC = 0.5;
  %INTERSECT_FRAC = 0.4;
  %MISS_FRAC = 0.75;
  %SUPP_PROB = 0.01;  
  
  %======== END SET PARAMETERS =======

  %======== STEP 2: NORMALIZE AFFINITIES AND GENERATE SIDE INFO =======
  [sA,D,sL,sM] = normalizeAffty(A);
  sqrtD    = sparse(diag(D.^0.5));
  sqrtDinv = sparse(diag(D .^ -0.5));
  % stationary distribution 
  u0 = D/sum(D);
  % DIFFUSE THE MARKOV MATRIX
  sMp = sM;      
  sMd = [];
  for k = 1:logpow 
    sMp = sMp * sMp;
    if k == (logpow-1)
      sMd = sMp;
    end
    if (logpow == 1)
      sMd = sM;
    end
  end
  %sLp = sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
  %======== END NORMALIZE AFFINITIES AND GENERATE SIDE INFO ===
  
  
  %======== STEP 3: KERNEL SELECTION ===========

  okIds = 1:size(sMp,2);
  boundIdsFlag = 0;
  if boundIdsFlag 
    % ids that are not on the boundary
    okIds = boundIds(sizeIm,1);
  end
  %okIds = 1:256;
  
  % a row vector
  mxM = MAX_HT_FRAC*max(sMp); 
  % same size as sMp: markov matrix diffused
  %gtM = sMp > (ones(size(sMp,1),1)*mxM); % aha time sink
  
  
  % saving memory with a loop
  %[sMpi,sMpj,sMps] = find(sMp);
  %sMps = sMps > mxM(sMpj);
  %gtM = sparse(sMpi,sMpj,sMps);
  
  
  %[mm,okIds] = sort(-mxM);
  
  % make a decision based on the stationary distribution
  [ss,okIds] = sort(-u0);
  
  % first pass: tight control on intersection fraction
  selectId = [];
  pixLoc   = zeros(prod(sizeIm),1);
  %pixLoc   = spalloc(prod(sizeIm),1,prod(sizeIm));  
  done = 0;
  k = 1;
  %nb = zeros(prod(sizeIm),1);
  %nb = spalloc(prod(sizeIm),1,prod(sizeIm));  
  while k <= length(okIds) & ~done
    
    % pick a kernel
    id = okIds(k);
    
    % pick nhbrs. these are atleast MAX_HT_FRAC*max(sMp).
    %nb = find(gtM(:,id));

    % doing the comparison where necessary    
    %nb = find(sMp(:,id) > mxM(id)); 
    nb = sMp(:,id) > mxM(id);     

    % what is common
    co = pixLoc .* (nb);
    %co = pixLoc .* sparse(nb);
    %co = sum(nb) - sum((nb - pixLoc) > 0) ;
    
    %if (length(co) < INTERSECT_FRAC*length(nb))
    %  selectId = [selectId id];
    %  pixLoc   = [pixLoc nb'];
    %end
    if (sum(co) < INTERSECT_FRAC*sum(nb))
    %if (co < INTERSECT_FRAC*sum(nb))      
      selectId = [selectId id];
      %pixLoc(find(nb)) = 1;
      pixLoc = pixLoc + nb;
    end
    
    k = k + 1;
  end
  %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
  %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
  pixMiss = find(1 - pixLoc);
  
  % second pass: relaxed control on intersection fraction  
  [mm,tt] = sort(-full(mxM(pixMiss)));
  % but based on stationary distribution
  %[mm,tt] = sort(-full(u0(pixMiss)));
  pixMiss = pixMiss(tt);
  
  skip = 0;
  if (~skip)
    k = 1;
    while k <= length(pixMiss) 
      % pick a kernel
      id = pixMiss(k);
      
      % pick nhbrs
      %nb = find(gtM(:,id));
      % doing the comparison where necessary	
      %nb = find(sMp(:,id) > mxM(id)); 
      %nb = sparse(sMp(:,id) > mxM(id));     
      nb = (sMp(:,id) > mxM(id));     
      
      % whats common
      co = pixLoc .* (nb);
      %co = pixLoc .* sparse(nb);      
      %co = sum(nb) - sum((nb - pixLoc) > 0) ;
      
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      if (sum(co) < MISS_FRAC*sum(nb))
      %if (co < INTERSECT_FRAC*sum(nb))      	
	selectId = [selectId id];
	%pixLoc(find(nb)) = 1;
	pixLoc = pixLoc + nb;
      end

      % OLD CODE
      % whats common
      %co = [];
      %co = intersect(pixLoc,nb);
      %if ( length(co) < MISS_FRAC*length(nb) )
      %	selectId = [selectId id];
      %	pixLoc   = [pixLoc nb'];
      %end
      
      k = k + 1;
    end
    %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
    %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
    pixMiss = find(1-pixLoc);
  end
  
  
  %% Diffusion Kernels, picked from sMd, not sMp
  %  where sMd is markov diffused to logpow-1.
  %K = sMd(:, selectId(:));
  % aug 10,'03 quick hack to see if sMp will work as well
  % in fact, it does seem to work. so i will leave it this way.
  K = sMd(:, selectId(:));
  % sort the kernels based on the pixel locations
  [selectId,id] = sort(selectId);
  % and shuffle the kernels
  K = K(:,id);
  %======== END KERNEL SELECTION =======

  %======== STEP 4: LATENT SPACE PROBABILITIES ======
  %[Ar,R,st,W,rbinNhbr] = latentProbs(K,u0);
  %[Ar,R,st,W,K] = setLatent(K,u0,1) ; % normalizes Ar by the median

  % set latent space. in particular, using kernels
  % K and the fine scale stationary distribution u0,
  % generate:
  %
  %  - responsibility matrix W 
  %  - latent space markov matrix R 
  %  - latent space stationary distribution st
  %  - latent space affinity matrix Ar.
  %
  % set small values of the markov matrix to zero using 
  % R + R', before generating the affinity matrix.
  % also scale the affinity matrix before returning.
  %
  % calls: emKernelFitNew
  %

  %%== OLD WAY
  %% W ownership matrix. last argument is display flag.
  %[st,W] = emKernelFit(u0,K,0); 
  %% latent space markov transition
  %R = W'*K; 
  %figure; plot(sum(R));
  %%== OLD WAY

  %%== NEW WAY
  [st,W,K,R,Ar] = emKernelFitNew(u0,K,0); % last value is display flag
  %%== NEW WAY

  ok = 0;
  
  if (ok)
  
    % thresholding small probabilities
    % consider: R = [a b; c d]
    % thresholding must be symmetric.  p(1->2) = b, is set to 0
    % only if p(2->1) = c, can also be set to 0.
    % hence use R + R' to come up with a Tolerance.
    % check if both p(1->2) and p(2->1) can be set to 0.
    flag  = 1;
    if (flag)
      trs = 0.001;
      % original line
      %Tol = 0.01*sum(R+R',2) * ones(1,size(R,2));
      % modified line
      Tol = trs*sum(R+R',2) * ones(1,size(R,2));
      RR  = (R + R') > 2*Tol;
      % find locations where p(1->2) is 0 but not p(2->1)
      RR(find(abs(RR-RR'))) = 1;
      R = R .* RR;
      Rs = sum(R,1);
      R = R ./ (ones(size(R,1),1)*Rs) ;
      % as R is quantized, update st, the latent space
      % stationary distribution by power iteration
      for k = 1:50
	st = R*st;
      end
      %figure; plot(sum(R)); 
      %fprintf('sum(st): %f \n',sum(st));
      
      % similarly, update the ownership matrix
      W = K * diag(st);
      W = diag(sum(W,2).^-1) * W;
    end
    
    % latent space affty
    %Ar = R * diag(st);
    Ar = R .* (st*ones(1,size(R,2)));
    
    % make sure it is symmetric 
    %Ar = (Ar + Ar')/2;
    
  else
    
    %Ar = R .* (st*ones(1,size(R,2)));

    trs = SUPP_PROB_FRAC;%0.01;
    %Tol = trs*sum(R+R',2) * ones(1,size(R,2));
    %RR  = (R + R') > 2*Tol;
    Tol = trs*sum(R+R',2);
    RR = (spdiags(Tol.^-1, 0, length(Tol), length(Tol)) * (R+R')) > 2;
    % find locations where p(1->2) is 0 but not p(2->1)
    RR(find(abs(RR-RR'))) = 1;
    
    Ar = Ar .* RR ;
    Dr = full(sum(Ar,1))';
    %R  = Ar .* (ones(size(Ar,1),1)*(Dr .^ -1));
    R  = Ar * spdiags(Dr .^ -1, 0, length(Dr), length(Dr));
    st = Dr/ sum(Dr);
    st = st(:);
    
    % similarly, update the ownership matrix
    W = K * spdiags(st, 0, length(st), length(st));
    W = spdiags(full(sum(W,2).^-1),0, size(W,1), size(W,1)) * W;
    
  end % check ok
  
  
  % scale latent space affinities. this will not affect
  % the transition matrix.
  Ar = Ar/median(full(sum(Ar,1)));

  rbinNhbr = Ar > 0;
  rbinNhbr = rbinNhbr - ...
      spdiags(ones(size(rbinNhbr,1),1),0,size(rbinNhbr,1),size(rbinNhbr,1));
  
  %======== END LATENT SPACE PROBABILITIES ======  
   
  return;
