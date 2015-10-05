function [st,own,K,R,Ar]  = emKernelFitNew_watershed(u0,K,L,p)
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
  %tol = 0.001;
  tol = 0.01;
  %TOT_ITR = 50; % ORIGINAL
  TOT_ITR = 50;
  
  if (nargin < 3)
    flag = 0;
  end

  % stationary distribution
  %u0 = D/sum(D);
  
  %st = ones(size(K,2),1); st = st/sum(st);
  st = u0'*K; st = st(:); st = st/sum(st);
  btsorig = KLdist(full(u0),full(K*st));

 

  %Q = zeros(size(K));
  for itr=1:TOT_ITR
    
    % E-step
    
    % copy this down below, too, once it works!!!
    own = K * spdiags(st,0,length(st),length(st));
    
    ownz = find(sum(own,2) == 0);
    
    % wait!! make sure that all rows in W sum to a nonzero. 
    pL = u0(p)/max(u0(p)); % a weight for the kernels
    % get pL on the same scale as own
    diff = max(pL)/max(own(:));
    pL = pL/diff;
    
    for ko = 1:size(ownz,1)
        k = ownz(ko);
        if mod(ko,100) == 1
            disp([num2str(ko),' : ',num2str(length(ownz))]);
        end
        if L(k) > 0
            w = L(k); % major kernel for this voxel
        else
            w = L(k-1);
        end
        we = pL(w);  % weight in this kernel
        nbweight = (1 - we); % divide remaining weight onto 6 neighbors
        nneigh = length(max(w-3,1):min(w+3,size(K,2))) - 1;
        ownx = k*ones(nneigh + 1,1); (k,max(w-3,1):min(w+3,size(K,2))) = nbweight/nneigh;
        owny = [max(w-3,1):min(w+3,size(K,2))];
        owns = nbweigt/nneigh;
    end
    
    Nown = sparse(ownx,owny,owns,size(own,1),size(own,2),length(owns));
    figure; spy(Nown);
    own = own + Nown;
    
    
    %own = diag(sum(own,2).^-1) * own;
    invsumown = sum(own,2).^-1 ;

    %own = own .* (invsumown*ones(1,size(own,2)));
    % this stmt is taking time. how about a for loop?
    %for zz = 1:size(own,2)
    %  own(:,zz) = own(:,zz) .* invsumown;
    %end
    own = spdiags(invsumown,0, size(own,1),size(own,1)) * own;
    
    %% own(i,k) is now the ownership of kernel k for data at pixel i

    % M-step
    st1 = own' * u0;

    % this works. 
    %K = diag(u0)*own*diag((st1+eps).^-1) ;
    % this is to speed up the computation 
    st2 = (st1+eps).^-1;
    
    % save time, rewrite this piece.
    %for zz = 1:size(own,2)
    %  K(:,zz) = u0 .* (own(:,zz)*st2(zz));
    %end
    
  % commented out 11-14
%     K = spdiags(u0, 0, size(own,1), size(own,1)) * ...
%         (own * spdiags(st2, 0, size(own,2), size(own,2)));
%     
    % will this translate to convergence based
    % on KL distance between u0 and K*st?
    step = max(abs(st1-st));
    st = st1;
    %pause;
    
    btsupd = KLdist(full(u0),full(K*st)) ;
    changeFrac = abs(btsupd-btsorig)/(btsorig+eps) ;
    

    if  changeFrac < tol
      break;
    end
    btsorig = btsupd ;
    
  end

  own = K * spdiags(st,0,length(st),length(st));
  %own = diag(sum(own,2).^-1) * own;
  invsumown = sum(own,2).^-1 ;
  % speed this up with a for loop !!
  if 1
    own = spdiags(invsumown,0, size(own,1),size(own,1)) * own;
    %for zz = 1:size(own,2)
    % own(:,zz) = own(:,zz) .* invsumown;
    %end
  end
   
  %% own(i,k) is now the ownership of kernel k for data at pixel i

  colones = ones(size(K,1),1) ;
  rowones = ones(1,size(K,2)) ;
  % get w(h|x)
  % speed this up with a for loop !!
  %Q = K.*(colones*st');
  Q = spalloc(size(K,1), size(K,2), size(K,2)*max(full(sum(K>0,1))));
  %for zz = 1:size(K,2)
  %  Q(:,zz) = K(:,zz)*st(zz);
  %end
  Q = K * spdiags(st,0,length(st), length(st));

  % save time. rewrite this.
  %Q = Q ./ (u0*rowones);
  % save more time. rewrite this as a for loop
  %Q = Q .* ((u0.^(-1))*rowones);  
  invu0 = u0.^-1;
  %for zz = 1:size(Q,2)
  %  Q(:,zz) = Q(:,zz) .* invu0;
  %end
  Q = spdiags(invu0, 0, length(invu0), length(invu0)) * Q;
  

  %for zz = 1:size(Q,2)
  %  Q(:,zz) = Q(:,zz) ./ u0;
  %end
  
  Q = Q';
  
  % get p(h'|h)
  %sumlatent = (Q .* (ones(size(Q,1),1)*u0'))*Q';
  %R = sumlatent ./ (ones(size(Q,1),1) * sum(sumlatent,1));
  Ar = (Q * spdiags(u0,0,length(u0),length(u0)))*Q';
  
  % data is really symmetric at this stage.
  % check out the suppression of small values in 
  % buildLatent.m or buildLatentProtein.m
  % after the call to emKernelFitNew.m
  %figure; spy(Ar-Ar'); max(max(abs(Ar-Ar')))/max(max(Ar))
  
  R  = Ar * spdiags(full(sum(Ar,1))'.^-1, 0, size(Ar,2), size(Ar,2));

  return;