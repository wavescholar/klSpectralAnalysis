function [aU,aS,Ur,Sr,itr] = coarseFineTPI(sL,Ar,K,nSvd,Ur,Sr,st)
%function [aU,aS] = coarseFine(sL,Ar,K, nSVD)
%function [aU,aS] = coarseFine(sL,Ar,K, nSVD,Ur,Sr,st)
%
% using latent space affinities Ar, interpolation 
% kernels K and fine scale normalized affinities sL, 
% generate fine scale eigenvectors/eigenvalues.
% use power iteration on the fine scale eigenvectors.
% use only the top nSvd eigenvectors. the eigenvalues
% are computed by the rayleigh function.
%
%% calls: gramFixed.m
% edited by vmb, 09-28-2010 --> changed mx

  USE_KIRCHHOFF = 0;
  USE_SPARSE_SVD = 0;  % 0 => use full svd on small problems, it is faster.
  POW_ITR_TR =  1e-4;%1e-3; % protein data: 1e-4
  TOL =  1e-4;%1e-3;
  TOTAL_POWER_ITR = 50;%50;%25;%25 and 20 work.
 
  %keyboard;
  %======== STEP 1:  ===========
  if nargin<5 | (prod(size(Ur)) == 0)

    % scaling Ar. not necessary, done by setLatent.m
    %Ar = Ar/median((sum(Ar,1)));
    
    Dr = sum(Ar,2);
    
    sqrtD    = spdiags(Dr.^0.5, 0, length(Dr), length(Dr));
    sqrtDinv = spdiags(Dr.^ -0.5, 0, length(Dr), length(Dr));
    
    if USE_KIRCHHOFF
      sLc = spdiags(Dr,0,length(Dr),length(Dr)) - Ar;
    else
      sLc = sqrtDinv*Ar*sqrtDinv;
    end
    
    if issparse(Ar) & USE_SPARSE_SVD
      [Ur,Sr,Vr] = svds(sLc,min(nSvd,size(Ar,2)));
      Sr = diag(Sr);
    else
      [Ur,Sr,Vr] = svd(full(sLc));
      Sr = diag(Sr);
    if 0  figure; plot(Sr);
    end % cancelled this plot (vmb, 2-3-12)
        % Crop to the first nSvd eigenvectors
      if nSvd < size(Ur,2)
	if USE_KIRCHHOFF
	  %Ur = Ur(:,end-nSvd+1:end);
	  %Sr = Sr(end-nSvd+1:end);	  
	  % they get sorted in the code below
	  % so this step is not necessary.
	  Ur = Ur(:,end:-1:end-nSvd+1);
	  Sr = Sr(end:-1:end-nSvd+1);	 
	  %Sr(1:10)
	else
	  Ur = Ur(:,1:nSvd);
	  Sr = Sr(1:nSvd);
	end
      end
    end
  end

  %%% DEBUG MODE 
  %keyboard;
  
  %======== END STEP 1  ===========

  %======== STEP 2:  ===========  
  % power iteration to update eigenvectors.
  % orthogonalize with fixed gram procedure.
  %
  % eigen values are obtaineed by rayleigh coefficient
  % measure. can the rayleigh coefficients be ever negative?
  %
  % to do: stopping criterion for power iterations
  % calls: gramFixed (to orthogonalize power iterated basis)
  %
  
  % interpolate, that is extend Ur to fine scale
  aU = K*Ur;
  %isreal(aU)

  
  if USE_KIRCHHOFF
  %  mx = max(Sr)
   % mx = 0.5*mean(sum(abs(sL),2));

    % this is what I have been using so far
   mx = 0.8*mean(sum(abs(sL),2));   % commented out by vmb, 09-28-10 

  %  mx = mean(sum(abs(sL),2));     % uncommented by vmb, 09-28-10  - does
  %  really well, except for EVectors - but the eigenvectors are fake,
  %  anyway, I think
 %   mx = median(sum(abs(sL),2));      % does well, but not as well as mean (except eigenvectors are better)      
 %  mx = max(sum(abs(sL),2))/2; % much smaller error - does best!!! (except for eigenvectors)
    
  %  mx = max(sum(abs(sL),2));    
    
    ml = size(sL,1);
    %keyboard;
   sL = spdiags(mx*ones(ml,1),0,ml,ml) - sL;
   
 %  sL = gershgorin(sL);
  
  end
  
  % orthogonalize right away (Apr'03)
  %[id,aU] = gramFixed(aU,0);    
  %aS = diag(aU'*sL*aU);
  %[aS,id] = sort(-aS);
  %aS = -aS;
  %aU = aU(:,id);

  % if one of the input arguments is the stationary 
  % distribution at the fine scale

  %if exist('st', 'var')
  %  aU(:,1) = sqrt(st(:));
  %end
  
  DEBUG = 0;
  itr = 1;
  if DEBUG
    fprintf('itr: ');  
  end

  Snsvd = min(Sr);
  TPI = min(TOTAL_POWER_ITR, floor(log(nSvd * eps/TOL)/log(Snsvd)));
  if TPI == 0
    TPI = 5;
  end
  
  %TPI = TOTAL_POWER_ITR;
  
  converged = 0; 
  %while (~converged & itr < 100)
  while (~converged)
    aUold = aU;
    %showIm(aU);

    for k = 1:TPI
      aU = sL*aU;
    end
    %isreal(aU)
    
    % orthogonalize in space
    %showIm(aU);
    [id,aU] = gramFixed(aU,0);    
    %showIm(aU);
    %fprintf('norm aU: %d %f %f\n',itr,max(sum(aU.*aU,1)),min(sum(aU.*aU,1)));
    aS = aU'*sL*aU;
    % orthogonalize in space residuals 
    [us,ss,vs] = svd(aS);
    % update in space eigenvectors
    aU = aU*us;
    %isreal(aU)
    % grab in space eigenvalues
    aS = diag(ss);
    
    [aS,id] = sort(-aS);
    aS = -aS;
    aU = aU(:,id);

    ok = 0;
    if ok
      % reset values in aS that are > 1
      aS(aS > 1) = 1.0;
      % set lengths of eigenvectors to be 1
      lengthU = sum(aU.*aU,1);
      lidx    = lengthU > 1.0;
      if (lidx)
	aU(:,lidx) = aU(:,lidx) ./ (ones(size(aU,1),1)*sqrt(lengthU(lidx))) ;
      end
    end
    
    innerP = 1 - sum(aU .* aU,1)';
    %innerP = 1 - abs(sum(aU .* aUold,1))';    
    %fprintf('innerP aU: %f %f\n',max(sum(aU.*aUold,1)),min(sum(aU.*aUold,1)));
    converged = max(abs(innerP)) < POW_ITR_TR;
    
    if DEBUG 
      pause;
    end

    itr = itr + 1;
    
  end
  %fprintf('itr: %f \n',itr);
  
  %======== END STEP 2  ===========  

  %fprintf('Coarse-Fine Iterations: %d\n',itr);
  
  
  return;
  

  %JUNK CODE

  %%ANGLE_TR = pi/36;
  %% the pow_itr_tr is roughly equal to (1 - cos(angle_tr))

  %GRAMS_ITR  = 50;
  %HALF_LIFE_TR = 50;

  
