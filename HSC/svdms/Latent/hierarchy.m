function [aU,aS] = hierarchy(A,sizeIm)

  debug = 0;

  A = sparse(A);
  D = sum(A,1)';                 % Normalize column sum to one.
  sM = A * sparse(diag(D.^-1));  % Markov
  sqrtD    = sparse(diag(D.^0.5));
  sqrtDinv = sparse(diag(D .^ -0.5));
  sL = sqrtDinv * A * sqrtDinv;

  logpow = 2; % 6 for the face image

  itr = 1;
  if (itr > 0)
    % pick kernels for latent space.
    latentPar.logpow   = logpow;
    latentPar.regular  = 1 ; % regular/irregular (1/0) sampling
    latentPar.minRes   = 0.005;
    latentPar.boundIdsFlag = 0; % ignore pixels on the boundaries
    latentPar.sizeIm   = sizeIm;
    [K,selectId] = pickLatent(sM,sqrtD,sqrtDinv,latentPar);
    % sort the kernels based on the pixel locations
    [selectId,id] = sort(selectId);
    % and shuffle the kernels
    K = K(:,id);
  else
    % kernel locations do not change. so diffuse fine scale 
    % markov matrix to pick out the kernels 
    sMd = sM;      
    for k = 1:logpow-1
      sMd = sMd * sMd;
    end
    sMp = sMd;
    %% Diffusion Kernels
    K = sMp(:,selectId);
  end  
  
  % set up latent space probabilities (stationary + markov)
  u0 = D/sum(D) ;
  % setLatent normalizes Ar by the median
  [Ar,R,st,W] = setLatent(K, u0,1) ; 
  rbinNhbr = Ar > 0;
  rbinNhbr = rbinNhbr - eye(size(rbinNhbr)) ;
  
  %% COARSE TO FINE INTERPOLATION
  [aU,aS,Ur,Sr] = coarseFine(sL,Ar,K) ;
  

  
  
  if (debug)

    % kernels with ids drawn as text
    kCols = ceil(selectId/sizeIm(1));
    kCols = kCols(:);
    kRows = selectId(:) - (kCols - 1)*sizeIm(1);
    kText = [];
    for t = 1:length(selectId)
      kText = strvcat(kText,sprintf('%d',t));
    end
    figure(114); clf;
    showIm(im); hold on;
    text(kCols,kRows,kText,'color','r','fontsize',12,'fontweight','demi');
    title(sprintf('Iteration: %d',itr),'fontsize',14);
    
    [Uf,Sf,Vf] = svd(full(sL));
    Sf = diag(Sf);
    figure(112);clf;
    plot(Sf(1:length(aS)),'bx-');hold on;
    plot(aS,'ro-');
    
    figure(115); clf;
    ct = 1;
    id2 = 6;
    for k = 1:id2
      subplot(2,id2,ct);
      showIm(reshape(aU(:,k),sizeIm));
      title(sprintf('Itp %d (%2.3f)',k,aS(k)));
      ct = ct + 1;
    end
    for k = 1:id2
      subplot(2,id2,ct);
      sgn = sign(aU(:,k)' * Uf(:,k));
      showIm(reshape(sgn*Uf(:,k),sizeIm));
      title(sprintf('%d (%2.3f)',k,Sf(k)));
      ct = ct + 1;
    end
  
  end
  
