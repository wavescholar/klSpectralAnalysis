function [aU,aS,Ur,Sr,itr] = coarseFine(sL,Ar,K,Ur,Sr,st)
%function [aU,aS] = coarseFine(sL,Ar,K)
%
% using latent space affinities Ar, interpolation 
% kernels K and fine scale normalized affinities sL, 
% generate fine scale eigenvectors/eigenvalues.
% use power iteration on the fine scale eigenvectors.
% use only the top 50 eigenvectors. the eigenvalues
% are computed by the rayleigh function.
%
%% calls: gramFixed.m
%

  POW_ITR_TR =  0.0001;
  ANGLE_TR = pi/36;
  % the pow_itr_tr is roughly equal to (1 - cos(angle_tr))
  
  GRAMS_ITR  = 50;
  HALF_LIFE_TR = 50;
  TOTAL_POWER_ITR = 25;
  
  %======== STEP 1:  ===========
  if ~exist('Ur', 'var') | (prod(size(Ur)) == 0)

    % scaling Ar. not necessary, done by setLatent.m
    %Ar = Ar/median((sum(Ar,1)));
    
    Dr = sum(Ar,2);
    
    sqrtD    = sparse(diag(Dr.^0.5));
    sqrtDinv = sparse(diag(Dr.^ -0.5));
    
    %[Ur,Sr,Vr] = svd( diag(Dr.^-0.5) * Ar * diag(Dr.^-0.5) );
    [Ur,Sr,Vr] = svd( full(sqrtDinv * Ar * sqrtDinv));
    
    Sr = diag(Sr);
    %Sr(1:10)
  end
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

  % orthogonalize right away (Apr'03)
  %[id,aU] = gramFixed(aU,0);    
  %aS = diag(aU'*sL*aU);
  %[aS,id] = sort(-aS);
  %aS = -aS;
  %aU = aU(:,id);

  % if one of the input arguments is the stationary 
  % distribution at the fine scale
  if exist('st', 'var')
    aU(:,1) = sqrt(st(:));
  end
  
  
  DEBUG = 0;
  itr = 1;
  if DEBUG
    fprintf('itr: ');  
  end

  converged = 0; 
  while (~converged)

    aUold = aU;
    
    for k = 1:TOTAL_POWER_ITR
      aU = sL*aU;
    end
    
    % orthogonalize 
    [id,aU] = gramFixed(aU,0);    
    aS = aU'*sL*aU;
    [us,ss,vs] = svd(aS);
    aU = aU*us;
    aS = diag(ss);
    
    [aS,id] = sort(-aS);
    aS = -aS;
    aU = aU(:,id);

    innerP = 1 - sum(aU .* aU,1)';
    converged = max(innerP) < POW_ITR_TR;
    
    if DEBUG 
      pause;
    end

    itr = itr + 1;
    
  end

  %======== END STEP 2  ===========  

  fprintf('Coarse-Fine Iterations: %d\n',itr);
  
  
  return;
  
  % junk code

    done = 0; % quick way to avoid going thro' this loop
  while (~done)

    aU = sL*aUold;

    aS = diag(aU'*sL*aU);
    [aS,id] = sort(-aS);
    aS = -aS;
    aU = aU(:,id);
    
    % begin angle test
    % cosine of angle
    cang = sum(aU.*aUold);

    nrm = sqrt(sum(aU.*aU));
    cang = cang ./ nrm;
    nrm = sqrt(sum(aUold.*aUold));
    cang = cang ./ nrm;

    if DEBUG
      figure(213);
      plot(cang,'x-');hold on;
    end

    % cos(0) = 1, cos(pi/2) = 0
    % vector that rotated the most will have
    % the miniumum abs(cang) value
    mang = min(abs(cang));
    if DEBUG
      fprintf('(%d %1.2f) ',itr,mang);
    end
    
    itr = itr + 1;
    
    if (abs(1-mang) < POW_ITR_TR)      
    %if (abs(mang) > cos(ANGLE_TR))
      done = 1;
    end 
    % end angle test on the eigenvectors

    
    % orthogonalize every once in a while
    if (rem(itr,GRAMS_ITR)==0) % 10 original value
      [id,aU] = gramFixed(aU,0);    
      aS = diag(aU'*sL*aU);
      [aS,id] = sort(-aS);
      aS = -aS;
      aU = aU(:,id);
    end

    if DEBUG 
      pause;
    end

    aUold = aU;
    
  end

  
  %[id,aU] = gramFixed(aU,1); % to debug
  [id,aU] = gramFixed(aU,0);
  aS = diag(aU'*sL*aU);
  [aS,id] = sort(-aS);
  aS = -aS;
  aU = aU(:,id);
