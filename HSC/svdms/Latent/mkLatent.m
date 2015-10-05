function [Ar,aU,aS] = mkLatent(sM,sL,D,sqrtD,sqrtDinv,sizeIm,logpow)
%function Ar = mkLatent(A,M,L,logpow)
% set up latent space  

  % for A, generate the corresponding
  % values for D, L and M
  %[sA,D,sL,sM] = normalizeAffty(A);
  %sqrtD = sparse(diag(D.^0.5));
  %sqrtDinv = sparse(diag(D .^ -0.5));
  
  %[U S V] = svds(sL,50); 
  %S = diag(S);
  %S(1:min(length(S),10))
  
  % generate latent space  
  regular = 0;
  if ~regular
    %% kernel selection: irregular pyramid
    %logpow = 2; % 4;
    sMp = sM;      
    for k = 1:logpow 
      sMp = sMp * sMp;
    end
    sLp =  sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
    minRes2 = 0.005;
    %% removing boundary effects while doing kernel selection.
    %% what are the ids along the four walls of the image?
    ids = boundIds(sizeIm,1);
    %% debug option on
    [selectId,c] = gramSelect(sLp(:,ids),minRes2,1);
    selectId = ids(selectId);
  else
    %% kernel selection: regular pyramid
    ok = 1;
    
    if (ok)
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
    else
      %% (2) sampling
      qq = reshape(1:prod(sizeIm),sizeIm);
      selectId = qq(1:2:size(qq,1),1:2:size(qq,2));
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
  
  %% mixing coefficient
  u0 = D/sum(D);
  
  % W ownership matrix
  [st,W] = emKernelFit(u0,K,0); 
  % latent space markov transition
  R = W'*K; 
  %figure; plot(sum(R));
  
  % latent space affty
  Ar = R * diag(st);
  % scaling of latent space affinities.
  % not necessary for this expt. however...
  Ar = Ar/mean((sum(Ar,1)));
  
  Dr = sum(Ar,2);
  [Ur,Sr,Vr] = svd( diag(Dr.^-0.5 ) * Ar * diag(Dr.^-0.5) );
  Sr = diag(Sr);
  %Sr(1:10)
  
  % interpolate, that is extend Ur to fine scale
  aU = K*Ur;
  
  % power iteration with gram selection
  %tt = 20;
  tt = min(50,size(aU,2));
  aU = aU(:,1:tt);

  % powIter calls gramFixed which has the debug option
  [aU,aS] = powIter(sL,aU);
  
  %for k = 1:2*tt
  %  aU = sL*aU;
  %end
  %[id,aU] = gramFixed(aU,1);
  %aS = diag(aU'*sL*aU);
  %[aS,id] = sort(-aS);
  %aS = -aS;
  %aU = aU(:,id);

