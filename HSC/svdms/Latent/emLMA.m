function [st,K,R]  = emLMA(jA,u0,K,st,flag)
%function [st,own,K,R]  = emLMA(jA,u0,K,st,flag)
% jA : joint distribution
% u0 : stationary distribution. model it as:  K*st
% K  : kernels for the mixture model
% st : mixture coefficients (stationary distribution
%      in the latent space)
% R  : transition matrix in the latent space
%
%% relating stationary distribution in the latent space st
%% to the stationary distribution in the full resolution
%% space u0
%% EM expression: u0' * log(K * st) + lambda * (sum(st) - 1)
  
  % very poor initialization
  %K = rand(size(K));
  %sumK = sum(K);
  %K = K./ (ones(size(K,1),1)*sumK) ;
  
  %tol = 0.0001;
  %tol = 0.001;
  tol = 0.01;
  
  if (nargin < 5)
    flag = 0;
  end
  
  % stationary distribution
  %u0 = D/sum(D);
  
  %st = ones(size(K,2),1); st = st/sum(st);
  %st = u0'*K; st = st(:); st = st/sum(st);
  st = st(:);
  colones = ones(size(K,1),1) ;
  rowones = ones(1,size(K,2)) ;
  st1     = zeros(length(st),1);
  K1      = sparse(zeros(size(K)));
  R       = sparse(zeros(size(K,2)));
  btsorig = KLdist(full(jA),full(K*diag(st)*K'));
  
  for itr = 0:50
    fprintf('%d...',itr);
    % E-step
    % p(h|x,x') = (g(x|h)g(x'|h)p(h))/(sum over h of the numerator)
    %   g(x|h)  = K 
    %   p(h)    = st
    %   p(x)    = u0
    % the numerator is a series of outer products and hence
    % lots of storage. so compute/store just the denominator
    % which is the normalizing factor for ownerships.
    sumown = (K.*(colones*st')) * K';
    nonzeroentries = find(sumown > 0);
    
    % M-step
    % slow
    slow = 1;
    if slow
      for idx = 1:length(st)
	tmp = K(:,idx)*K(:,idx)'*st(idx);
	tmp(nonzeroentries) = tmp(nonzeroentries) ./ sumown(nonzeroentries); 
	tmp = tmp .* jA;
	st1(idx) = sum(sum(tmp)) ;
	K1(:,idx) = sum(tmp,2)/st1(idx);
      end
    else
      % fast 
      for idx = 1:length(st)
	% how fast can i can run the loop above
      end
    end    

    % slow
    %btsupd = KLdist(full(jA),full(K1*diag(st1)*K1'));
    % fast ?
    btsupd = KLdist(full(jA),full( (K1 .* (colones*st1')) * K1'));
    changeFrac = abs(btsupd-btsorig)/(btsorig+eps) ;
    
    if (flag)
      figure(112); clf; 
      subplot(2,1,1)
      plot(st,'go-'); 
      hold on
      stem(st1,'b'); 
      subplot(2,1,2);
      plot(u0,'b-');    
      hold on;
      plot(K1*st1,'go-');
      title(sprintf('itr: %d KL: %2.3f KLchange: %2.3f',itr,btsupd,changeFrac));
    end
    st = st1;
    K  = K1;
    

    if  changeFrac < tol
      break;
    end
    btsorig = btsupd ;
    
    %pause;
  end

  % get w(h|x)
  Q = K .* (colones*st');
  Q = Q ./ (u0*rowones);
  Q = Q';
  
  % get R = p(h'|h)
  % slow version
  %sumlatent = Q*diag(u0)*Q';
  % fast version?
  sumlatent = (Q.*(ones(size(Q,1),1)*u0'))*Q';
  R         = sumlatent ./ ( ones(size(Q,1),1) * sum(sumlatent,1) );
  
  fprintf('\n');
