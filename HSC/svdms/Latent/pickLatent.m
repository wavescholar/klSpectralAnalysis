function [K,selectId,sMp] = pickLatent(sM,sqrtD,sqrtDinv,latentPar,st)
%function [K,selectId] = pickLatent(sM,sqrtD,sqrtDinv,latentPar)
% diffuse Markov matrix to 2^logpow iterations
% and pick kernels for the latent space.
% i avoid pixels along the boundary. not exactly sure why.
%
% right now this code is called with regular = 1
% and option 3 is hardcoded.
    
  logpow  = latentPar.logpow;
  regular = latentPar.regular;
  minRes2 = latentPar.minRes;
  sizeIm  = latentPar.sizeIm;
  boundIdsFlag = latentPar.boundIdsFlag;

  % DIFFUSE THE MARKOV MATRIX
  %logpow = 2; % 4;
  sMp = sM;      
  for k = 1:logpow 
    sMp = sMp * sMp;
  end
  sLp = sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
  
  % generate latent space  
  if ~regular
    %% kernel selection: irregular pyramid
    %minRes2 = 0.005;

    %% removing boundary effects while doing kernel selection.
    %% what are the ids along the four walls of the image?
    ids = 1:size(sLp,2);
    if boundIdsFlag 
      ids = boundIds(sizeIm,1);
    end
    
    %% debug option on/off
    debug = 0;
    [selectId,c] = gramSelect(sLp(:,ids),minRes2,debug);
    selectId = ids(selectId);
    
  else
    
    %% kernel selection: regular pyramid
    ok = 3;
    
    if (ok==1)
      %% sqrt(2) sampling
      selectId = [];
      step = 2;
      for k = 0:step:sizeIm(2)-step+1
	selectId  = [selectId ([1:2:sizeIm(1)] + k*sizeIm(1))];
      end
      for k = 1:step:sizeIm(2)-step+1
	selectId  = [selectId ([2:2:sizeIm(1)] + k*sizeIm(1))];
      end
      selectId = sort(selectId);
    end

    if (ok ==2)
      %% (2) sampling
      qq = reshape(1:prod(sizeIm),sizeIm);
      selectId = qq(1:2:size(qq,1),1:2:size(qq,2));
      selectId = sort(selectId(:));
    end

    if (ok == 3)

      %okIds = 1:size(sLp,2);
      okIds = 1:size(sMp,2);
      boundIdsFlag = 0;
      if boundIdsFlag 
	% ids that are not on the boundary
	okIds = boundIds(sizeIm,1);
      end
      %okIds = 1:256;

      % find an ordering on the kernels based on max height.
      % increasing this # will increase the total 
      % number of kernels selected
      maxHtFrac = 0.5;%0.6;%0.5
      % length of the intersection 
      intersectFrac = 0.4;%0.5;%0.51;%0.5;%
      missFrac = 0.85;%0.75;%0.85

      % a row vector
      mxM = maxHtFrac*max(sMp); 
      % same size as sMp: markov matrix diffused
      gtM = sMp > (ones(size(sMp,1),1)*mxM); 
      [mm,okIds] = sort(-mxM);

      [ss,okIds] = sort(-st);
      
      selectId = [];
      pixLoc   = [];
      done = 0;
      k = 1;
      while k <= length(okIds) & ~done

	% pick a kernel
	id = okIds(k);
	
	% pick nhbrs
	nb = find(gtM(:,id));
	% whats common
	co = [];
	co = intersect(pixLoc,nb);
	if (length(co) < intersectFrac*length(nb))
	  selectId = [selectId id];
	  pixLoc   = [pixLoc nb'];
	end
	k = k + 1;
      end
      pixLoc  = intersect(1:prod(sizeIm),pixLoc);
      pixMiss = setdiff(1:prod(sizeIm),pixLoc);

      [mm,tt] = sort(-mxM(pixMiss));
      pixMiss = pixMiss(tt);
      k = 1;
      while k <= length(pixMiss) 
	% pick a kernel
	id = pixMiss(k);
	% pick nhbrs
	nb = find(gtM(:,id));
	% whats common
	co = [];
	co = intersect(pixLoc,nb);
	if ( length(co) < missFrac*length(nb) )
	  selectId = [selectId id];
	  pixLoc   = [pixLoc nb'];
	end
	k = k + 1;
      end
      pixLoc  = intersect(1:prod(sizeIm),pixLoc);
      pixMiss = setdiff(1:prod(sizeIm),pixLoc);
      
      %selectId = sort([selectId pixMiss]);

    end
    
    %pp = zeros(sizeIm);
    %pp(selectId) = 1;
    %figure; showIm(pp)
    
  end

  % diffuse fine scale markov matrix
  % to pick out the kernels 
  sMd = sM;      
  %logpow = 2;
  %for k = 1:2*logpow-1
  for k = 1:logpow-1
    sMd = sMd * sMd;
  end
  sMp = sMd;
  
  %% Diffusion Kernels
  K = sMp(:, selectId(:));
  

