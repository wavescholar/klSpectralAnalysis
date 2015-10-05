function [Ar,R,st,W,K] = setLatent(K,u0,flag)
%function [Ar,R,st,W] = setLatent(K,u0)
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
% calls: emKernelFit
%

  % flag to do thresholding of markov matrix
  if nargin < 3
    flag = 0;
  end
  
  %%== OLD WAY
  %% W ownership matrix. last argument is display flag.
  %[st,W] = emKernelFit(u0,K,0); 
  %% latent space markov transition
  %R = W'*K; 
  %figure; plot(sum(R));
  %%== OLD WAY
  

  %%== NEW WAY
  [st,W,K,R,Ar] = emKernelFitNew(u0,K,0); 
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
      trs = 0.1;
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
    Ar = (Ar + Ar')/2;
    
  else
    
    %Ar = R .* (st*ones(1,size(R,2)));

    trs = 0.01;
    Tol = trs*sum(R+R',2) * ones(1,size(R,2));
    RR  = (R + R') > 2*Tol;
    % find locations where p(1->2) is 0 but not p(2->1)
    RR(find(abs(RR-RR'))) = 1;
    
    Ar = Ar .* RR ;
    Dr = sum(Ar,1);
    R  = Ar .* (ones(size(Ar,1),1)*(Dr .^ -1));
    st = Dr/ sum(Dr);
    st = st(:);
    
    % similarly, update the ownership matrix
    W = K * diag(st);
    W = diag(sum(W,2).^-1) * W;
    
  end % check ok
  
  
  % scale latent space affinities. this will not affect
  % the transition matrix.
  Ar = Ar/median((sum(Ar,1)));

  return;
  
  
  
  old = 0;
  if (old)
    [st,W,K,R] = emKernelFitNew(u0,K,0); 
    
    % thresholding small probabilities
    % consider: R = [a b; c d]
    % thresholding must be symmetric.  p(1->2) = b, is set to 0
    % only if p(2->1) = c, can also be set to 0.
    % hence use R + R' to come up with a Tolerance.
    % check if both p(1->2) and p(2->1) can be set to 0.
    if (flag)
      Tol = 0.01*sum(R+R',2) * ones(1,size(R,2));
      RR  = (R + R') > 2*Tol;
      % find locations where p(1->2) is 0 but not p(2->1)
      RR(find(abs(RR-RR'))) = 1;
      R = R .* RR;
      Rs = sum(R,1);
      R = R ./ (ones(size(R,1),1)*Rs) ;
      % as R is quantized, update st, the latent space
      % stationary distribution by power iteration
      for k = 1:20
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
    Ar = (Ar + Ar')/2;
    
  end
  
  
