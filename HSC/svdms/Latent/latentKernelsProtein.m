function [selectId,newprotInfo] = latentKernelsProtein(D,sMp,protInfo)

  selectId = [];
  newprotInfo = protInfo;
  
  % work on each sub-unit
  for i = 1:protInfo.subunits
  
    % which chains belong to the sub-unit i
    chainIds = find(protInfo.chains == i);

    % find the first chain and its boundaries
    firstChain = min(chainIds);
    if firstChain == 1
      loC = 1;
    else
      loC = protInfo.chLen(firstChain-1)+1;
    end
    hiC = protInfo.chLen(firstChain);

    % TO DO: diffuse the markov matrix for the firstChain
    
    % kernel selection on the firstChain. call to localKernels
    % code is below.
    localId = localKernels(D(loC:hiC),sMp(loC:hiC,loC:hiC));
    localId = find(localId);
    
    % using symmetry pick the suitable ids on all chains
    % of the sub-unit
    for j = 1:length(chainIds)
      if chainIds(j) == 1
	loC = 0;
      else
	loC = protInfo.chLen(chainIds(j)-1);
      end
      selectId = [selectId; loC + localId(:)];
      newprotInfo.chLen(chainIds(j)) = length(localId);
    end
  end
  %keyboard;
  newprotInfo.chLen = cumsum(newprotInfo.chLen);
  
  % TO DO: diffuse the markov matrix for the selectId columns
  
  return;
  
%================================================================
function [selectId] = localKernels(D,sMp)  
%=== kernel selection parameters ===
% find an ordering on the kernels based on max height.
% increasing this # will increase the total 
% number of kernels selected
    
  %Aug10, '04 for protein data  
  % original values: [0.55, 0.4 , 0.75]
    
  % BEST VALUES (JUNE '05): [0.5 0.3 0.75] 
  MAX_HT_FRAC = 0.5;%0.5; % 
  % length of the intersection 
  INTERSECT_FRAC = 0.3;%0.3;
  MISS_FRAC = 0.75;
  
  sizeIm = size(sMp,1);
  %======== STEP 3: KERNEL SELECTION ===========
  okIds = 1:size(sMp,2);
  
  % a row vector 
  %mxM = MAX_HT_FRAC*max(sMp); 
  % BIG CHANGE HERE. punt for now.
  mxM = MAX_HT_FRAC*(diag(sMp)); 
  mxM = mxM';
  
  % same size as sMp: markov matrix diffused
  %gtM = sMp > (ones(size(sMp,1),1)*mxM); % aha time sink
  
  
  % saving memory with a loop
  %[sMpi,sMpj,sMps] = find(sMp);
  %sMps = sMps > mxM(sMpj);
  %gtM = sparse(sMpi,sMpj,sMps);
  
  
  %[mm,okIds] = sort(-mxM);
  
  % make a decision based on the stationary distribution
  %[ss,okIds] = sort(-D);
  % BIG CHANGE HERE
  [ss,okIds] = sort(-diag(sMp));  
  
  % first pass: tight control on intersection fraction
  %selectId = [];
  selectId   = 1:prod(sizeIm) < 0;
  selectId   = selectId(:);
  
  %pixLoc   = zeros(prod(sizeIm),1);
  pixLoc   = 1:prod(sizeIm) < 0;
  pixLoc   = pixLoc(:);
  
  done = 0;
  k = 1;
  %nb = zeros(prod(sizeIm),1);
  nb = 1:prod(sizeIm) < 0;
  kernelCount = 1;
  while k <= length(okIds) & ~done
    
    % pick a kernel
    id = okIds(k);
    %if id == 8
    %  keyboard;
    %end
    
    if ~pixLoc(id)
      % pick nhbrs. these are atleast MAX_HT_FRAC*max(sMp).
      %nb = find(gtM(:,id));
      
      % doing the comparison where necessary    
      %nb = find(sMp(:,id) > mxM(id)); 
      nb = sMp(:,id) >= mxM(id);     
      
      % whats common
      %co = pixLoc .* nb;
      co = pixLoc & nb;
      
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      
      % sum(co) gives what is common with the already
      % selected nodes. sum(nb) gives the numer of
      % nodes that pass the threshold test. 
      % node id is included if what is common
      % is less than the numer of nodes that pass
      % the threshold test.
      if (sum(co) < INTERSECT_FRAC*sum(nb))
	%selectId = [selectId id];
	
	% original code
	selectId(id) = 1;
	% new code
	%selectId(id) = kernelCount;      
	%kernelCount = kernelCount + 1;
	
	pixLoc((nb)) = 1;
      end
      
    end
    k = k + 1;
  end
  %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
  %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
  pixMiss = find(1 - pixLoc);
  %length(pixMiss)
  
  %keyboard;
  
  % second pass: relaxed control on intersection fraction  
  [mm,tt] = sort(-full(mxM(pixMiss)));
  % but based on stationary distribution
  %[mm,tt] = sort(-full(u0(pixMiss)));
  %[mm,tt] = sort(-full(D(pixMiss)));  % because u0 is not defined here.
  %pixMiss = pixMiss(tt);

  % skip = 0. checking if i should skip the following.
  skip = 0;
  if (~skip & length(pixMiss > 0))
    k = 1;
    while k <= length(pixMiss) 
      % pick a kernel
      id = pixMiss(k);
      
      % pick nhbrs
      %nb = find(gtM(:,id));
      % doing the comparison where necessary	
      %nb = find(sMp(:,id) > mxM(id)); 
      nb = sMp(:,id) > mxM(id);     
      
      % whats common
      %co = pixLoc .* nb;
      co = pixLoc & nb;
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      if (sum(co) < MISS_FRAC*sum(nb))
	%selectId = [selectId id];

	% original code
	selectId(id) = 1;
	% new code
	%selectId(id) = kernelCount;      
	%kernelCount = kernelCount + 1;
	
	pixLoc((nb)) = 1;
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
  %length(pixMiss)  
  
  %======== END KERNEL SELECTION =======

%================================================================
%% ORIGINAL CODE. TRYING SOME NEW THINGS IN THE CODE ABOVE
%================================================================
function [selectId] = localKernelsOriginal(D,sMp)  
%=== kernel selection parameters ===
% find an ordering on the kernels based on max height.
% increasing this # will increase the total 
% number of kernels selected
    
  %Aug10, '04 for protein data  
  % original values: [0.55, 0.4 , 0.75]
    
  % BEST VALUES (JUNE '05): [0.5 0.3 0.75] 
  MAX_HT_FRAC = 0.5;%0.5; % 
  % length of the intersection 
  INTERSECT_FRAC = 0.3;%0.3;
  MISS_FRAC = 0.75;
  
  sizeIm = size(sMp,1);
  %======== STEP 3: KERNEL SELECTION ===========
  okIds = 1:size(sMp,2);
  
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
  [ss,okIds] = sort(-D);
  
  % first pass: tight control on intersection fraction
  %selectId = [];
  selectId   = 1:prod(sizeIm) < 0;
  selectId   = selectId(:);
  
  %pixLoc   = zeros(prod(sizeIm),1);
  pixLoc   = 1:prod(sizeIm) < 0;
  pixLoc   = pixLoc(:);
  
  done = 0;
  k = 1;
  %nb = zeros(prod(sizeIm),1);
  nb = 1:prod(sizeIm) < 0;
  kernelCount = 1;
  while k <= length(okIds) & ~done
    
    % pick a kernel
    id = okIds(k);
    if ~pixLoc(id)
      % pick nhbrs. these are atleast MAX_HT_FRAC*max(sMp).
      %nb = find(gtM(:,id));
      
      % doing the comparison where necessary    
      %nb = find(sMp(:,id) > mxM(id)); 
      nb = sMp(:,id) >= mxM(id);     
      
      % whats common
      %co = pixLoc .* nb;
      co = pixLoc & nb;
      
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      
      if (sum(co) < INTERSECT_FRAC*sum(nb))
	%selectId = [selectId id];
	
	% original code
	selectId(id) = 1;
	% new code
	%selectId(id) = kernelCount;      
	%kernelCount = kernelCount + 1;
	
	pixLoc((nb)) = 1;
      end
      
    end
    k = k + 1;
  end
  %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
  %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
  pixMiss = find(1 - pixLoc);
  length(pixMiss)
  
  % second pass: relaxed control on intersection fraction  
  [mm,tt] = sort(-full(mxM(pixMiss)));
  % but based on stationary distribution
  %[mm,tt] = sort(-full(u0(pixMiss)));
  %[mm,tt] = sort(-full(D(pixMiss)));  % because u0 is not defined here.
  %pixMiss = pixMiss(tt);

  % skip = 0. checking if i should skip the following.
  skip = 0;
  if (~skip & length(pixMiss > 0))
    k = 1;
    while k <= length(pixMiss) 
      % pick a kernel
      id = pixMiss(k);
      
      % pick nhbrs
      %nb = find(gtM(:,id));
      % doing the comparison where necessary	
      %nb = find(sMp(:,id) > mxM(id)); 
      nb = sMp(:,id) > mxM(id);     
      
      % whats common
      %co = pixLoc .* nb;
      co = pixLoc & nb;
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      if (sum(co) < MISS_FRAC*sum(nb))
	%selectId = [selectId id];

	% original code
	selectId(id) = 1;
	% new code
	%selectId(id) = kernelCount;      
	%kernelCount = kernelCount + 1;
	
	pixLoc((nb)) = 1;
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
  length(pixMiss)  
  
  %======== END KERNEL SELECTION =======
