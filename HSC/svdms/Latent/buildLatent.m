function [sL,Ar,K,R,st,W,rbinNhbr,selectId,sMd,newW] = buildLatent(A,logpow,sizeIm,fctr,oldW,kfctr,lev,smallprob) 
%
% buildLatent: set up latent space, build kernels, 
%              generate transitional matrix for 
%              the low-dim space.
%
% calls to:
%   normalizeAffty.m 
%   emKernelFitNew.m
%
% fctr ... added by vmb, 11-12-10: controls lower limit for kernels to
% select (smallest kernel must have stationary prob >= max(kernel stationary prob)*fctr
% oldW ... added by vmb, 11-17-10: for suppressing kernels: each kernel
% must contain a small fraction of the finest data (to suppress kernels
% that have practically no members)


  %======== STEP 1: SET PARAMETERS ===========

  %=== kernel selection parameters ===
  % find an ordering on the kernels based on max height.
  % increasing this # will increase the total 
  % number of kernels selected
  MAX_HT_FRAC = .55;%0.25;%0.55; % 0.55 as of Aug10, '04 for protein data
  % length of the intersection 
  INTERSECT_FRAC = 0.1;%0.4; % lower numbers create unsmooth divisions, and odd separations
  MISS_FRAC = 0.85; %.75
  % to threshold small probabilities. look in STEP 4 below.
  SUPP_PROB_FRAC = smallprob; % 0.01;
  
  %% parameter set that works: Apr 23, '03
  %MAX_HT_FRAC = 0.5;
  %INTERSECT_FRAC = 0.4;
  %MISS_FRAC = 0.75;
  %SUPP_PROB = 0.01;  
  
  %======== END SET PARAMETERS =======

  
  %======== STEP 2a: NORMALIZE AFFINITIES AND GENERATE SIDE INFO =======
  [sA,D,sL,sM,sH] = normalizeAffty(A); % sH added by virginia, 09-24 
  
  
  minlev = 300;
  if ((lev >= minlev) && (size(sL,1) < 3000))
      
      %=========== STEP 2b: cut sensitive edges ==================
      % added by vmb, 11-18-2010
     
      tA = sA; % working matrix for each iteration
      figure; spy(tA);
      
      %%% GET LARGEST CONNECTED COMPONENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % make sure you have just one connected component
      
      % choose a threshold
      nz = find(tA>0);
      mintA = min(tA(nz));
      maxtA = max(tA(nz));
      breadthtA = maxtA - mintA;
      tauf = .05 + (lev-minlev)*.01; % favorite: either .05 or .01...
      tau = mintA + tauf*breadthtA;
      
      %tau = 0;
      Y = gingcca_speedup(tA,tau); %cca(sparse(tA),tau,1);
      z = max(Y(:));
      
%       % take a look at components:
%       lbl_Y = gingcca_speedup(tA,tau);
%       imagesc(lbl_Y);
      
      disp([num2str(z),' components found']);
      %       if z > size(tA,1)*.25
      %           disp('move on from here: the matrix is already very segmented');
      %
      %       else
      
      comps = cell(z,1);
      csizes = zeros(z,1);
      for cn = 1:z
          comps{cn} = find(Y == cn);
          csizes(cn) = size(comps{cn},1);
      end
      % don't actually need to sort!!! could save time..
      % [ac bc] = sort(csizes,'descend');
      % comps = comps(bc);
      % csizes = ac;
      
      % check if the largest component is large enough to consider individually
      % ummm... consider all components that make up at least a percent of
      % the data.
      % at higher levels, let a lower percentage be okay.
      %perc = .2 + (lev - minlev)*.05
      isolated = max(size(tA,1)/90,2);
      sigs = find(csizes >= isolated);% size(tA,1)*perc);
      nsig = size(sigs,1);   
      disp(['there are ',num2str(nsig),' significant components']);
      
      %% before moving on, remove (set to very low prob) components that are completely
      %% disconnected from cutA. 
      %isolated = size(tA,1)/90;
      disconnect = find(csizes < isolated);
      figure; spy(tA);
      title(['points ignored bc < ',num2str(isolated),' connected points']);
      rs = zeros(size(tA,1),1);
      ndisc = 0;
      for discomp = disconnect'
          c = comps{discomp};
         % cutA(c) = 1e-8; % a small number
          rs(ndisc+1:ndisc+length(c)) = c;
          ndisc = ndisc + length(c);

      end
      rs = rs(1:ndisc);
      [x y] = ind2sub([size(tA,1),size(tA,2)],rs);
      r = sparse(x,y,ones(ndisc,1),size(tA,1),size(tA,2),ndisc);
      hold on; spy(r,'r',10);
      
      
      %% Now analyze connected components! 
      % old goal: remove sensitive edges
      % New goal: choose kernels for each connected component!
      % then let each component contribute a percentage of its kernels,
      % based on its size, and the ranking of its kernels
      
     selectIds  = cell(nsig,1); %zeros(size(A,1),1);
     u0_selectIds = cell(nsig,1); 
     nIds = 0; 
      
     % rearrange indicing of connected components
     lbl = diag(Y);
     compsL = cell(z,1);
     csizesL = zeros(z,1);
     for cn = 1:z
         compsL{cn} = find(lbl == cn);
         csizesL(cn) = size(compsL{cn},1);
     end
     
     for si = 1:nsig
          sig = sigs(si);
          
          c = comps{sig};
          nc = size(c,1);
         
          % plot the connected component:
          figure; spy(tA);
          [x y] = ind2sub([size(tA,1),size(tA,2)],c);
          r = sparse(x,y,tA(c),size(tA,1),size(tA,2),nc);
          hold on;
          spy(r,'r',20); title('current component that we will choose kernels from');
          clear r; % don't actually need this
          
          c = compsL{sig};
          nc = size(c,1);
          
          nA = tA(c,c);
        %  figure; spy(nA); title('nA','FontSize',20);
          [ttA,ttD,ttL,ttM,ttH] = normalizeAffty(nA);
          sizeImt = [size(ttA,1) 1];
          
          % how many cc does this have?
          [~, z2] = gingcca(nA,tau); % [lbl z] = cca(sparse(tA),tau,1);
          disp([num2str(z2),' components found']);
          if z2 > 1
              disp('this should only have one component now...');
          end
          
          %%%%%%% Choose some kernels from this component! %%%%%%%%
          
          u0 = ttD/sum(ttD);
          % DIFFUSE THE MARKOV MATRIX   (8.3.1 in thesis)
          tMp = ttM;
          tMd = [];
          for k = 1:logpow
              tMp = tMp * tMp;
              if k == (logpow-1)
                  tMd = tMp;
              end
              if (logpow == 1)
                  tMd = ttM;
              end
          end
          
          %   %======== STEP 3: KERNEL SELECTION =========== (8.3.3 in thesis)
          
          okIds = 1:size(tMp,2);
          boundIdsFlag = 0;
          if boundIdsFlag
              % ids that are not on the boundary
              okIds = boundIds(sizeImt,1);
          end
          mxM = MAX_HT_FRAC*max(tMp);
          
          % make a decision based on the stationary distribution
          [~,okIds] = sort(-u0);
          selectIdt   = 1:prod(sizeImt) < 0;
          selectIdt   = selectIdt(:);
          pixLoc   = 1:prod(sizeImt) < 0;
          pixLoc   = pixLoc(:);
          
          done = 0;
          k = 1;
          nb = 1:prod(sizeImt) < 0;
          kernelCount = 1;
          while k <= length(okIds) & ~done
              % pick next kernel (kernels are sorted by size(u0), and removed when
              % assigned elsewhere
              id = okIds(k);
              nb = tMp(:,id) > mxM(id);    % all rows within .55*max(this column)
              co = pixLoc & nb;
              if (sum(co) < INTERSECT_FRAC*sum(nb))
                  selectIdt(id) = 1;
                  pixLoc((nb)) = 1;
              end
              k = k + 1;
          end
          pixMiss = find(1 - pixLoc);
          
          [~,tt] = sort(-full(mxM(pixMiss)));
          pixMiss = pixMiss(tt);
          
          skip = 0;
          if (~skip)
              k = 1;
              while k <= length(pixMiss)
                  % pick a kernel
                  id = pixMiss(k);
                  nb = tMp(:,id) > mxM(id);
                  co = pixLoc & nb;
                  if (sum(co) < MISS_FRAC*sum(nb))
                      selectIdt(id) = 1;
                      pixLoc((nb)) = 1;
                  end
                  k = k + 1;
              end
              pixMiss = find(1-pixLoc);
          end
          
          
          %% Diffusion Kernels, picked from sMd, not sMp
          
          K = tMd(:, selectIdt);     
          selectIdt = find(selectIdt);
          
          % sort selectId by u0
          u0selectIds = u0(selectIdt);
          [au bu] = sort(u0selectIds,'descend');
          u0selectIdts = au;
          selectIdt = selectIdt(bu); 
          
          u0_selectIds{si} = u0selectIdts;
          selectIds{si} = c(selectIdt); % convert back from this component to the full matrix!
          
          nIds = nIds + size(selectIdt,1);
          
          
          
%           %%% FIND SENSITIVE EDGES WITH EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%
%           
%           [tU tS] = svds(tL,10);
%           ts = diag(tS);
%           tsqrtDinv = tD.^-0.5; %spdiags(D.^.5,0,length(D),length(D));
%           nev = min(10,length(ts)); %number of eigenvectors to look at
%           
%           zev = find(ts == 0);
%           if isempty(zev) == 0 % so if there are EVs = 0
%               ts'
%               nev = min(find(ts == 0)) - 1;
%               tU = tU(:,1:nev);
%               ts = ts(1:nev);
%               disp('there are eigenvalues = 0. is this okay? only using non-zeros.');
%           end
%           
%           % plot eigenvectors!!
%           nsub = ceil(sqrt(nev));
%           figure;
%           for id = 1:nev
%               subplot(nsub,nsub,id);
%               plot(tU(:,id)); axis xy; axis tight;
%               title(['eigenvector ',num2str(id)]);
%           end
%           
%           ts = min(ts, 1-10*eps); % don't want log(< 1)
%           halfL = -log(2) * log(ts.^2).^-1 ;
%           beta0 = halfL(end);
%           
%           % record indices (x,y) and labels (idxs) of points in lower triangle of A.
%           % points with low class labels (max = 32) will be removed.
%           xs = cell(nev,1);
%           ys = cell(nev,1);
%           idxs = cell(nev,1);
%           for id = 2:nev
%               dLogHalf = perturbBeta(beta0,ts(id),tU(:,id),tsqrtDinv);
%               % fix the diagonal
%               dLogHalf = dLogHalf .* (tA>0);
%               [xs{id} ys{id} idxs{id}] = plotContactColor(dLogHalf,sprintf('%d',id));
%               %dLogHalfAll = [dLogHalfAll ; dLogHalf(find(tA>0))];
%           end
%           
%           % look at histogram assignments neatly
%           figure;
%           for id = 2:nev
%               subplot(3,3,id-1);
%               hist(idxs{id},32);
%           end
%           
%           %%% CUT EDGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           
%           % cut all edges with idx{id} < 20 (range from 1-32).
%           % do this for each EV looked at
%           % Make cuts in original matrix A!!! (including sections that were removed
%           % bc disconnected)
%                     
%           figure; spy(cutA);
%           r = zeros(size(cutA));
%           hold on;
%           
%           cutlim = 21;% + (lev-3)*1;
%           for id = 2:nev
%               cut = find(idxs{id} < cutlim);
%               % x,y,id
%               cutx = xs{id}(cut);
%               cuty = ys{id}(cut);
%               for t = 1:size(cut,1)
%                   cutA(c(cutx(t)),c(cuty(t))) = 0;
%                   cutA(c(cuty(t)),c(cutx(t))) = 0;
%                   
%                   % show where cuts are in map = this is just for
%                   % debugging!!
%                   r(c(cutx(t)),c(cuty(t))) = 1;
%                   r(c(cuty(t)),c(cutx(t))) = 1;
%               end
%              % r = sparse(c(cutx),c(cuty),ones(length(cutx)),size(nA,1),size(nA,2),length(cutx));
%               spy(r,'r',10);
%               
%           end
%           disp('paint components that were cut here, to see if I am doing this part right!');
          
     end
     
   %  figure; spy(cutA);
     
     %=========== STEP 2c: renormalize affinities ==================
     % added by vmb, 11-18-2010
   %  [sA,D,sL,sM,sH] = normalizeAffty(cutA); % sH added by virginia, 09-24
  
  % now let each component contribute its highest ranked kernels to selectId.
  % the number of kernels per component is chosen based on the component's percent
  % size.
  
  if nsig > 0
  
  psizes = csizes(sigs)/sum(csizes(sigs));
  kerfact = 1;%.9; % only keep 95% of total number of kernels found from all components
  nkernels = round(nIds*kerfact);
  nIds = 0;
  
  selectId = zeros(nkernels,1);
  for si = 1:nsig
      nksig = floor(psizes(si)*nkernels);
      % choose the nkisg highest ranked kernels (based on u0) to add to selectId.
      % make sure that sig has enough kernels to offer.
      sigkerns = length(selectIds{si});
      fromsig = selectIds{si}(1:min(nksig,sigkerns));
      if ~isempty(fromsig)
          lfromsig = size(fromsig,1);
          selectId(nIds+1:nIds+lfromsig) = fromsig;
          nIds = nIds + lfromsig;
      end
  end
  selectId = selectId(1:nIds);
  selectId = sort(selectId);
  
  
  else
      selectId = 1:size(D,1);
  end
  
  u0 = D/sum(D);
  % DIFFUSE THE MARKOV MATRIX   (8.3.1 in thesis)
  sMp = sM;
  sMd = [];
  for k = 1:logpow
      sMp = sMp * sMp;
      if k == (logpow-1)
          sMd = sMp;
      end
      if (logpow == 1)
          sMd = sM;
      end
  end
  
  
  K = sMd(:, selectId);    
  [st,W,K,R,Ar] = emKernelFitNew(u0,K,0); % doesn't update kernels!!!
  
%   % update kernels slightly:
%   %     1) each kernel leader must have a minimal probability
%         Kstatprob = u0(selectId);
%         figure; subplot(2,1,1);
%         [as bs] = sort(selectId);
%         plot(as,Kstatprob(bs)); % the stat. dist. of the kernel leaders
%         maxKstat = max(Kstatprob);
%         toosmall = selectId(find(Kstatprob < maxKstat*fctr));
%         bigenough = selectId(find(Kstatprob >= maxKstat*fctr));
%         hold on; plot(bigenough,zeros(length(bigenough),1),'xr')
%         plot(selectId,maxKstat*fctr,'m--');
%         
%         if size(toosmall,1) > 0
%           disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc stat dist too small']); 
%           selectId = bigenough;
%           K = sMd(:, bigenough);
%          [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);  
%          
%         end   
%   
%   %     2) each kernel must have at least two? members from finest layer
%     % use W to suppress kernels with low membership
%         newW = oldW*W;
%         [a c] = max(newW,[],1); 
%         subplot(2,1,2);
%         [as bs] = sort(selectId);
%         plot(as,a(bs));
%        % plot(selectId,a);
%         min(a)
%         bigenough = selectId(find(a > kfctr));
%         hold on; plot(bigenough,zeros(length(bigenough),1),'xr')
%         toosmall = selectId(find(a <= kfctr));
%         if size(toosmall,1) > 0
%             disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc highest % membership too small']); 
%             selectId = bigenough;
%             K = sMd(:, bigenough);
%             [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);  
%             newW = oldW*W;
%         end
  newW = [];
  
  else % end Steps 2b-c (only done for lev >= minlev)
      
      
      
      
      
      
      %===========================================================
      %% ANALYZING SHIFTED KIRCHOFF MATRIX. REWRITING THE ARRAY SL
      % estimate the upperbound on the max eigenvalue from Gershgorin
      % circles as in the paper by ACE. however, it does not give
      % a good approximation for the higher order vectors
      % here are few attempts at fixing it.
      %  mx = max(sum(abs(sH),2));
      % mx = 1.2*mean(sum(abs(sH),2));
      
      % mx = 0.5*mean(sum(abs(sH),2));  % works the best among these
      % sL = spdiags(mx*ones(length(D),1),0,length(D),length(D)) - sH;
      
      % from ACE (virginia, 09-28-10)
      % we want to get the smallest EVs of sL, and this method gives us the
      % largest. So we have to flip sL.
      
      % sL = spdiags(D,0,length(D),length(D)) - sA;
      
      % sL = gershgorin(sL);
      
      
      %===========================================================
      
      % changed to save memory by vmb, 10-6-10 - these variables aren't even
      % used, anyway, so hardly matter and take up lots of memory
      % and commented out by vmb, 11-8-10
      % sqrtD  = spdiags(D .^ 0.5,0,length(D),length(D)); %sparse(diag(D.^0.5));
      % sqrtDinv = spdiags(D .^ -0.5,0,length(D),length(D)); %sparse(diag(D .^ -0.5));
      % stationary distribution
      u0 = D/sum(D);
      % DIFFUSE THE MARKOV MATRIX   (8.3.1 in thesis)
      sMp = sM;
      sMd = [];
      for k = 1:logpow
          sMp = sMp * sMp;
          if k == (logpow-1)
              sMd = sMp;
          end
          if (logpow == 1)
              sMd = sM;
          end
      end
      %sLp = sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5
      %======== END NORMALIZE AFFINITIES AND GENERATE SIDE INFO ===
      
      
      %======== STEP 3: KERNEL SELECTION =========== (8.3.3 in thesis)
      
      okIds = 1:size(sMp,2);
      boundIdsFlag = 0;
      if boundIdsFlag
          % ids that are not on the boundary
          okIds = boundIds(sizeIm,1);
      end
      %okIds = 1:256;
      
      % a row vector
      mxM = MAX_HT_FRAC*max(sMp);
      % same size as sMp: markov matrix diffused
      %gtM = sMp > (ones(size(sMp,1),1)*mxM); % aha time sink
      
      
      % saving memory with a loop
      %[sMpi,sMpj,sMps] = find(sMp);
      %sMps = sMps > mxM(sMpj);
      %gtM = sparse(sMpi,sMpj,sMps);
      
      
      %[mm,okIds] = sort(-mxM);
      
      % make a decision based on the stationary distribution
      [~,okIds] = sort(-u0);
      
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
      while k <= length(okIds) & ~done
          
          % add: if pixLoc(k) == 0 do, otherwise continue (k = k+1) or something...
          
          % pick biggest okID as first kernel (kernels are sorted by size(u0), and removed when
          % assigned elsewhere
          id = okIds(k);
          
          % pick nhbrs. these are atleast MAX_HT_FRAC*max(sMp).
          %nb = find(gtM(:,id));
          
          % doing the comparison where necessary
          %nb = find(sMp(:,id) > mxM(id));
          nb = sMp(:,id) > mxM(id);    % all rows within .55*max(this column)
          
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
      %pixLoc  = intersect(1:prod(sizeIm),pixLoc);
      %pixMiss = setdiff(1:prod(sizeIm),pixLoc);
      pixMiss = find(1 - pixLoc);
      
      % this could be simplified, and the parameter could be entirely
      % removed!! (miss_frac = . * intersec_frac)...
      % this was to make sure that everybody was covered by a kernel
      % how often does this even matter? 
      disp(['pixMiss: ',num2str(length(pixMiss))]);
      
      % second pass: relaxed control on intersection fraction
      [mm,tt] = sort(-full(mxM(pixMiss)));
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
      
      % is anybody still missing???
      disp(['pixMiss after second round: ',num2str(length(pixMiss))]);
      
      
      %% Diffusion Kernels, picked from sMd, not sMp
      %  where sMd is markov diffused to logpow-1.
      %K = sMd(:, selectId(:));
      % aug 10,'03 quick hack to see if sMp will work as well
      % in fact, it does seem to work. so i will leave it this way.
      %K = sMd(:, selectId(:));
      K = sMd(:, selectId);
      
      
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
      [st,W,K,R,Ar] = emKernelFitNew(u0,K,0); % last value is display flag
      %%== NEW WAY
      
 
  
  % update kernels slightly:
  %     1) each kernel leader must have a minimal probability
        Kstatprob = u0(selectId);
        figure; subplot(3,1,1);
        plot(selectId,Kstatprob); % the stat. dist. of the kernel leaders
        maxKstat = max(Kstatprob);
        toosmall = selectId(find(Kstatprob < maxKstat*fctr));
        bigenough = selectId(find(Kstatprob >= maxKstat*fctr));
        hold on; plot(bigenough,zeros(length(bigenough),1),'xr')
        plot(selectId,maxKstat*fctr,'m--');
        
        if size(toosmall,1) > 0
          disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc stat dist too small']); 
          selectId = bigenough;
          K = sMd(:, bigenough);
         [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);  
         
        end   
  
  %     2) each kernel must have at least two? members from finest layer
    % use W to suppress kernels with low membership
        newW = oldW*W;
        [a c] = max(newW,[],1); 
        subplot(2,1,2);
        plot(selectId,a);
        min(a)
        bigenough = selectId(find(a > kfctr));
        hold on; plot(bigenough,zeros(length(bigenough),1),'xr')
        toosmall = selectId(find(a <= kfctr));
        if size(toosmall,1) > 0
            disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc highest % membership too small']); 
            selectId = bigenough;
            K = sMd(:, bigenough);
            [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);  
            newW = oldW*W;
        end
 %newW = []; 
 
 
    % Each KERNEL MUST HAVE SOME MEMBERS!!! (otherwise what's the point??)
%     newW = oldW*W;
%     snewW = size(newW,2);
%     [a cln] = max(newW,[],2);
%     
%     if snewW < 50
%         colors = rand(snewW,3);
%         subplot(3,1,2); hold on;
%         for j = 1:snewW
%             plot(newW(:,j),'Color',colors(j,:));
%         end
%     end
%     
%     toosmall = selectId(setdiff([1:snewW],unique(cln)));
%    
%     if size(toosmall,1) > 0
%         bigenough = selectId(intersect([1:snewW],unique(cln)));
%         disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc no members']);
%         selectId = bigenough;
%         subplot(3,1,3);
%         [a2 b2] = max(W,[],2);
%         plot(a2);
%         hold on; plot(toosmall,zeros(length(toosmall),1),'xr')
%         title('kernels removed bc no members');
%        
%         K = sMd(:, bigenough);
%         [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);
%         newW = oldW*W;
%         snewW = size(newW,2);
%         
%         if snewW < 50
%             colors = rand(snewW,3);
%             subplot(3,1,3); hold on;
%             for j = 1:snewW
%                 plot(newW(:,j),'Color',colors(j,:));
%             end
%             title('new kernels!!');
%         end
%     end
 
    
    
    
%     %% EXTREME VERSION: each kernel must have a certain percent (> .5) of
%     %% members
%       % Each KERNEL MUST HAVE SOME MEMBERS!!! (otherwise what's the point??)
%     kfctr = .6;
%     newW = oldW*W;
%     snewW = size(newW,2);
%     [a cln] = max(newW,[],2);
%     
%     if snewW < 50
%         colors = rand(snewW,3);
%         subplot(3,1,2); hold on;
%         for j = 1:snewW
%             plot(newW(:,j),'Color',colors(j,:));
%         end
%     end
%     
%     toosmall = selectId(find(max(newW,[],1) <= kfctr));
%    
%     if size(toosmall,2) > 0
%         bigenough = selectId(find(max(newW,[],1) > kfctr));
%         disp([num2str(size(toosmall,1)),' of ',num2str(size(selectId,1)),' kernels left out bc no members']);
%         selectId = bigenough;
%         subplot(3,1,3);
%         [a2 b2] = max(W,[],2);
%         plot(a2);
%         hold on; plot(toosmall,zeros(length(toosmall),1),'xr')
%         title(['kernels removed bc never has at least ',num2str(kfctr),'% membership']);
%        
%         K = sMd(:, bigenough);
%         [st,W,K,R,Ar] = emKernelFitNew(u0,K,0);
%         newW = oldW*W;
%         snewW = size(newW,2);
%         
%         if snewW < 50
%             colors = rand(snewW,3);
%             subplot(3,1,3); hold on;
%             for j = 1:snewW
%                 plot(newW(:,j),'Color',colors(j,:));
%             end
%             title('new kernels!!');
%         end
%     end
    
    

  end  % end if lev < minlev
        
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
    RR = (spdiags(Tol.^-1, 0, length(Tol), length(Tol)) * (R+R')) > 2;
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
    RR = sparse(iii,jjj,sss);

    Ar = Ar .* RR ;

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

  rbinNhbr = Ar > 0;
  rbinNhbr = rbinNhbr - ...
      spdiags(ones(size(rbinNhbr,1),1),0,size(rbinNhbr,1),size(rbinNhbr,1));
  
  %======== END LATENT SPACE PROBABILITIES ======  
   
  return;
