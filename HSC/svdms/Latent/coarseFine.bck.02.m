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

  POW_ITR_TR =  1e-6;%0.0001;
  TOTAL_POWER_ITR = 50;%25 and 20 work.
  
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
    
    % orthogonalize in space
    [id,aU] = gramFixed(aU,0);    
    aS = aU'*sL*aU;
    % orthogonalize in space residuals 
    [us,ss,vs] = svd(aS);
    % update in space eigenvectors
    aU = aU*us;
    % grab in space eigenvalues
    aS = diag(ss);
    
    [aS,id] = sort(-aS);
    aS = -aS;
    aU = aU(:,id);

    % reset values in aS that are > 1
    aS(aS > 1) = 1.0;
    
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
  

  %JUNK CODE

  %%ANGLE_TR = pi/36;
  %% the pow_itr_tr is roughly equal to (1 - cos(angle_tr))

  %GRAMS_ITR  = 50;
  %HALF_LIFE_TR = 50;

  
