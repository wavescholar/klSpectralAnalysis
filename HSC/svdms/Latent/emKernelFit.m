function [st,own]  = emKernelFit(u0,K,flag)
%function [st,own] = emKernelFit(u0,K,flag)
% u0 : stationary distribution. model it as:  K*st
% K  : kernels for the mixture model
% st : mixture coefficients (stationary distribution
%      in the latent space)
% own: ownership matrix
%
%% relating stationary distribution in the latent space st
%% to the stationary distribution in the full resolution
%% space u0
%% EM expression: u0' * log(K * st) + lambda * (sum(st) - 1)
  
  
  %tol = 0.0001;
  tol = 0.01;
  
  if (nargin < 3)
    flag = 0;
  end
  
  % stationary distribution
  %u0 = D/sum(D);
  
  %st = ones(size(K,2),1); st = st/sum(st);
  st = u0'*K; st = st(:); st = st/sum(st);
  btsorig = KLdist(full(u0),full(K*st));
  
  for itr=0:75
    own = K * diag(st);
    
    %own = diag(sum(own,2).^-1) * own;
    invsumown = sum(own,2).^-1 ;
    own = own .* (invsumown*ones(1,size(own,2)));
    
    %% own(i,k) is now the ownership of kernel k for data at pixel i
    st1 = own' * u0;
    
    % will this translate to convergence based
    % on KL distance between u0 and K*st?
    step = max(abs(st1-st));
    st = st1;
    %if  step < tol
    %  break;
    %end
    
    btsupd = KLdist(full(u0),full(K*st));
    changeFrac = abs(btsupd-btsorig)/btsorig;
    %% plot the convergence
    if (flag)
      figure(111); clf; 
      subplot(2,1,1)
      plot(st,'go-'); 
      hold on
      stem(st1,'b'); 
      subplot(2,1,2);
      plot(u0,'b-');    
      hold on;
      plot(K*st1,'go-');
      title(sprintf('itr: %d KL: %2.3f KLchange: %2.3f',itr,btsupd,changeFrac));
    end
    if  changeFrac < tol
      break;
    end
    btsorig = btsupd ;
    
    %pause;
  end

  own = K * diag(st);
  %own = diag(sum(own,2).^-1) * own;
  invsumown = sum(own,2).^-1 ;
  own = own .* (invsumown*ones(1,size(own,2)));
  
  %% own(i,k) is now the ownership of kernel k for data at pixel i
