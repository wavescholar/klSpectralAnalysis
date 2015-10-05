function [K,selectId,sMp] = latentKernels(sM,sqrtD,sqrtDinv,logpow,sizeIm)
%
% actual work done in: pickLatent.m
%
  
  %logpow = 2;
  itr = 1;
  if (itr > 0)
    % pick kernels for latent space.
    latentPar.logpow   = logpow;
    latentPar.regular  = 1 ; %regular/irregular (1/0) sampling
    latentPar.minRes   = 0.005;
    latentPar.boundIdsFlag = 0; % 1 = ignore pixels on the boundaries
    latentPar.sizeIm   = sizeIm;

    st = sqrtD.^2;
    st = st/sum(st);
    [K,selectId,sMp] = pickLatent(sM,sqrtD,sqrtDinv,latentPar,st);

    % sort the kernels based on the pixel locations
    [selectId,id] = sort(selectId);
    % and shuffle the kernels
    K = K(:,id);
  else
    % diffuse fine scale markov matrix
    % to pick out the kernels 
    sMd = sM;      
    for k = 1:logpow-1
      sMd = sMd * sMd;
    end
    sMp = sMd;
    %% Diffusion Kernels
    K = sMp(:,selectId);
  end  
