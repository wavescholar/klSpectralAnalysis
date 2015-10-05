function [sL,Ar,K,R,st,W,rbinNhbr,selectId,sMd] = buildLatent_orig_msm2(A) 
% this version uses sMp = sMd*sMd (requires more memory, but hopefully
% faster kernel selection than buildLatent_orig_msm.m)
%
% buildLatent: set up latent space, build kernels, 
%              generate transitional matrix for 
%              the low-dim space.
%
% calls to:
%   normalizeAffty.m 
%   emKernelFitNew.m
%
  %======== STEP 1: SET PARAMETERS ===========

  logpow = 2; 
  sizeIm = size(A,1); 
  
  %=== kernel selection parameters ===
  % find an ordering on the kernels based on max height.
  % increasing this # will increase the total 
  % number of kernels selected
  MAX_HT_FRAC = .55%.1%.35%0.55;%0.5; % 0.55 as of Aug10, '04 for protein data
  % length of the intersection 
  INTERSECT_FRAC = .4%.95%0.4;%0.4;
  MISS_FRAC = .75%.9%9%0.75;
  % to threshold small probabilities. look in STEP 4 below.
  SUPP_PROB_FRAC = 0.01;
  
  %% parameter set that works: Apr 23, '03
  %MAX_HT_FRAC = 0.5;
  %INTERSECT_FRAC = 0.4;
  %MISS_FRAC = 0.75;
  %SUPP_PROB = 0.01;  
  
  %======== END SET PARAMETERS =======

  %======== STEP 2: NORMALIZE AFFINITIES AND GENERATE SIDE INFO =======
  [sA,D,sL,sM] = normalizeAffty(A);
  disp('A normalized');
  
  %===========================================================
  %% ANALYZING SHIFTED KIRCHOFF MATRIX. REWRITING THE ARRAY SL
  % estimate the upperbound on the max eigenvalue from Gershgorin
  % circles as in the paper by ACE. however, it does not give
  % a good approximation for the higher order vectors
  % here are few attempts at fixing it. 
  %mx = max(sum(abs(sH),2)); 
  %mx = 1.2*mean(sum(abs(sH),2)); 

  %mx = 0.5*mean(sum(abs(sH),2));  % works the best among these
  %sL = spdiags(mx*ones(length(D),1),0,length(D),length(D)) - sH;
  
  % sL = spdiags(D,0,length(D),length(D)) - sA;
  
  %===========================================================
    
%  sqrtD    = sparse(diag(D.^0.5));
 % sqrtDinv = sparse(diag(D .^ -0.5));
  % stationary distribution 
  u0 = D/sum(D);
%   % DIFFUSE THE MARKOV MATRIX
%   sMp = sM;      
%   sMd = [];
%   for k = 1:logpow 
%       k
%     sMp = sMp * sMp;
%     if k == (logpow-1)
%       sMd = sMp;
%     end
%     if (logpow == 1)
%       sMd = sM;
%     end
%   end
  %sLp = sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
  
  % vmb, 6-13-2011. changed so that sMp = sMd*sMd is never fully computed
  disp('diffusing sM twice');
  tic;
  %sMd = sM*sM; % diffuse once here. diffuse a second time individually for each column while picking kernels
  sMp = mpower(sM,4); %sMp = mpower2(sM,4);
  toc;
  disp('sM diffused twice');
  
  
  %======== END NORMALIZE AFFINITIES AND GENERATE SIDE INFO ===
  
  
  %======== STEP 3: KERNEL SELECTION ===========

  okIds = 1:size(sMp,2);
  boundIdsFlag = 0;
  if boundIdsFlag 
    % ids that are not on the boundary
    okIds = boundIds(sizeIm,1);
  end
  %okIds = 1:256;
  
  % a row vector
  mxM = MAX_HT_FRAC*max(sMp); 
% vmb, 6-13-2011. changed so that sMp = sMd*sMd is never fully computed
  
  % same size as sMp: markov matrix diffused
  %gtM = sMp > (ones(size(sMp,1),1)*mxM); % aha time sink
  
  
  % saving memory with a loop
  %[sMpi,sMpj,sMps] = find(sMp);
  %sMps = sMps > mxM(sMpj);
  %gtM = sparse(sMpi,sMpj,sMps);
  
  
  %[mm,okIds] = sort(-mxM);
  
  % make a decision based on the stationary distribution
  [ss,okIds] = sort(-u0);
  
  % first pass: tight control on intersection fraction
  %selectId = [];
  selectId   = 1:prod(sizeIm) < 0;
  selectId   = selectId(:);
  
  %pixLoc   = zeros(prod(sizeIm),1);
  pixLoc   = 1:prod(sizeIm) < 0;
  pixLoc   = pixLoc(:);
  
  done = 0;
  k = 1;
  %nb = zeros(prod(sizeIm),1);
  nb = 1:prod(sizeIm) < 0;
  kernelCount = 1;
  disp('beginning kernel selection');
  while k <= length(okIds) & ~done
      if mod(k,500) == 1
          disp([num2str(k),': ',num2str(length(okIds))]);
      end
    
    % pick a kernel
    id = okIds(k);
    
    % pick nhbrs. these are atleast MAX_HT_FRAC*max(sMp).
    %nb = find(gtM(:,id));

    % doing the comparison where necessary    
    %nb = find(sMp(:,id) > mxM(id)); 
    
    % vmb, 6-13-2011. changed so that sMp = sMd*sMd is never fully computed
  %  sMp = sMd*sMd(:,id);
  %  nb = sMp > MAX_HT_FRAC*max(sMp); %mxM(id);    
    nb = sMp(:,id) > mxM(id);     


    % whats common
    %co = pixLoc .* nb;
    co = pixLoc & nb;
    
    %if (length(co) < INTERSECT_FRAC*length(nb))
    %  selectId = [selectId id];
    %  pixLoc   = [pixLoc nb'];
    %end
    if (sum(co) < INTERSECT_FRAC*sum(nb))
      %selectId = [selectId id];
      
      % original code
      selectId(id) = 1;
      % new code
      %selectId(id) = kernelCount;      
      %kernelCount = kernelCount + 1;

      pixLoc((nb)) = 1;
    end
    
    k = k + 1;
  end
  
  disp('finished first round of kernel selection');
  %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
  %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
  pixMiss = find(1 - pixLoc);
  disp(['pixMiss: ',num2str(length(pixMiss))]);
  
  % second pass: relaxed control on intersection fraction  
  [mm,tt] = sort(-full(mxM(pixMiss))); 
%  mxD = MAX_HT_FRAC*max(sMd);  % use this instead, bc didn't compute mxM (bc no sMp)
%  [mm,tt] = sort(-full(mxD(pixMiss))); clear mxD;
  
  % but based on stationary distribution
  %[mm,tt] = sort(-full(u0(pixMiss)));
  pixMiss = pixMiss(tt);
  
  skip = 0;
  if (~skip)
    k = 1;
    while k <= length(pixMiss) 
      % pick a kernel
      id = pixMiss(k);
      
      % pick nhbrs
      %nb = find(gtM(:,id));
      % doing the comparison where necessary	
      %nb = find(sMp(:,id) > mxM(id)); 
      nb = sMp(:,id) > mxM(id);     
      
    % vmb, 6-13-2011. changed so that sMp = sMd*sMd is never fully computed
   % sMp = sMd*sMd(:,id);
   % nb = sMp > MAX_HT_FRAC*max(sMp); %mxM(id);    
      
      % whats common
      %co = pixLoc .* nb;
      co = pixLoc & nb;
      %if (length(co) < INTERSECT_FRAC*length(nb))
      %  selectId = [selectId id];
      %  pixLoc   = [pixLoc nb'];
      %end
      if (sum(co) < MISS_FRAC*sum(nb))
          %selectId = [selectId id];
          
          % original code
          selectId(id) = 1;
          % new code
          %selectId(id) = kernelCount;
          %kernelCount = kernelCount + 1;
          
          pixLoc((nb)) = 1;
      end

      % OLD CODE
      % whats common
      %co = [];
      %co = intersect(pixLoc,nb);
      %if ( length(co) < MISS_FRAC*length(nb) )
      %	selectId = [selectId id];
      %	pixLoc   = [pixLoc nb'];
      %end
      
      k = k + 1;
    end
    %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
    %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
    pixMiss = find(1-pixLoc);
  end
  
   disp(['pixMiss: ',num2str(length(pixMiss))]);
  
  %% Diffusion Kernels, picked from sMd, not sMp
  %  where sMd is markov diffused to logpow-1.
  %K = sMd(:, selectId(:));
  % aug 10,'03 quick hack to see if sMp will work as well
  % in fact, it does seem to work. so i will leave it this way.
  %K = sMd(:, selectId(:));
  clear sMp;
  disp('diffusing sM once');
  tic;
  sMd = mpower(sM,2);%= mpower2(sM,2);
  toc; 
  K = sMd(:, selectId);  
  
  % figure(82); clf;  plot(sum(K,1)); title([' sum K across columns. min row sum = ',num2str(min(sum(K,2)))],'FontSize',18); 
 disp('computed reduced diffused markov matrix, K = sMd*sMd(:,selectId)');

  % AVOID SORTING TO RETAIN THE KERNEL SELECTION ORDER
  % sort the kernels based on the pixel locations
  %[selectId,id] = sort(selectId); % clearly wrong
  % UNCOMMENT THE NEXT TWO OPERATIONS IF SORTED IS NEEDED
  %[selectId,id] = sort(find(selectId));  
  selectId = find(selectId);
  %qqq = find(selectId); [qqs,qqi] = sort(-u0(qqq));selectId=qqq(qqi);
  % and shuffle the kernels
  %K = K(:,id);
  %size(K)
  
  %======== END KERNEL SELECTION =======

  %======== STEP 4: LATENT SPACE PROBABILITIES ======
  %[Ar,R,st,W,rbinNhbr] = latentProbs(K,u0);
  %[Ar,R,st,W,K] = setLatent(K,u0,1) ; % normalizes Ar by the median

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
  % calls: emKernelFitNew
  %

  %%== OLD WAY
  %% W ownership matrix. last argument is display flag.
  %[st,W] = emKernelFit(u0,K,0); 
  %% latent space markov transition
  %R = W'*K; 
  %figure; plot(sum(R));
  %%== OLD WAY

  %%== NEW WAY
  disp('running em');
  [st,W,K,R,Ar] = emKernelFitNew_msm(u0,K,0); % last value is display flag
   disp('finished em');
   
 % trying to fix memory problems, 06 - 18- 2011
  clear W;
  u0 = 0;
 sizAr = size(Ar,2);
 if sizAr < 20000
     R  = Ar * spdiags((sum(Ar,1))'.^-1, 0, sizAr, sizAr); % R is the markov transition matrix..
 else
     % threshold!
     whos Ar
     [iii jjj sss] = find(Ar > max(max(Ar))/3000);
     Ar = sparse(iii,jjj,sss,size(Ar,1),size(Ar,2),length(sss));
     whos Ar
     ssA =sum(Ar,2);
     length(find(ssA == 0)), 
     clear ssA;
   %  Ar = Ar(find(Ar > max(max(Ar))/100));
     R  = Ar * spdiags((sum(Ar,1))'.^-1, 0, sizAr, sizAr); % R is the markov transition matrix..
 end
  
  disp('finished emKernelFitNew');
  disp(['size R: ',num2str(size(R,1)),' x ',num2str(size(R,2))]);
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
              trs = 0.001;
              % original line
              %Tol = 0.01*sum(R+R',2) * ones(1,size(R,2));
              % modified line
              Tol = trs*sum(R+R',2) * ones(1,size(R,2));
              RR  = (R + R') > 2*Tol;
              size(RR)
              % find locations where p(1->2) is 0 but not p(2->1)
              length(find(abs(RR-RR')))
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
          %Ar = (Ar + Ar')/2;
          
      else
          
          %Ar = R .* (st*ones(1,size(R,2)));
          
          trs = SUPP_PROB_FRAC;%0.01;
          %Tol = trs*sum(R+R',2) * ones(1,size(R,2));
          %RR  = (R + R') > 2*Tol;
          Tol = trs*sum(R+R',2);
  
          disp('computing RRn, RR');
          RRn = spdiags(Tol.^-1, 0, length(Tol), length(Tol)) * (R+R');
          RR = RRn > 2;
         % RR = (spdiags(Tol.^-1, 0, length(Tol), length(Tol)) * RRp) > 2;
          disp('computed RR part 1');
          
          % added by vmb, 06-15-2011, to avoid rows in A/R/W that contain only 0.
          RRnrowmax = max(RRn,[],2);
          clear RRn;
          % fix zero rows of RR
          zrow = find(RRnrowmax < 2);
          nzrow = size(zrow,1);
          RR = RR + sparse(zrow,zrow,ones(nzrow,1),size(RR,1),size(RR,2),nzrow);
          clear RRnrowmax; 
          clear zrow nzrow;
          disp('added ones to diag of zero rows in RR');
          
          
          % find locations where p(1->2) is 0 but not p(2->1)
          %size(RR)
          [rri,rrj,rrs] = find(abs(RR-RR'));
          %[length(xyz) max(xyz)]
          % this should do, but it stopped working when N = 512
          %RR(find(abs(RR-RR'))) = 1;
          [iii,jjj,sss] = find(RR);
          
          iii = [rri ; iii];
          jjj = [rrj ; jjj];
          sss = [ones(length(rrs),1) ; sss];
          %[min(iii) min(jjj) min(sss)]
          RR = sparse(iii,jjj,sss,size(RR,1),size(RR,2),size(iii,1)); % specifying dimensions in RR fixes bug in following multiplication
          
          disp('multiplying');
          Ar = Ar .* RR ;
          disp('finished multiplying');
          
          %% HACK ALERT. June '05 for protein data
          % make sure it is symmetric
          Ar = (Ar + Ar')/2;
          
          Dr = full(sum(Ar,1))';
          %R  = Ar .* (ones(size(Ar,1),1)*(Dr .^ -1));
          R  = Ar * spdiags(Dr .^ -1, 0, length(Dr), length(Dr));
          st = Dr/ sum(Dr);
          st = st(:);
          
          % similarly, update the ownership matrix
          W = K * spdiags(st, 0, length(st), length(st));
          W = spdiags(full(sum(W,2).^-1),0, size(W,1), size(W,1)) * W;
          
      end % check ok
      
  % scale latent space affinities. this will not affect
  % the transition matrix.
  Ar = Ar/median(full(sum(Ar,1)));

%  rbinNhbr = Ar > 0;
%  rbinNhbr = rbinNhbr - ...
%      spdiags(ones(size(rbinNhbr,1),1),0,size(rbinNhbr,1),size(rbinNhbr,1));
rbinNhbr = []; % vmb, 6-14-2011 bc of memory issues with Ar > 0
  
  %======== END LATENT SPACE PROBABILITIES ======  
   
  return;
