%==============================================
% instead of going thro' the code here
% i wrote: exptsLatent.m    (Feb 03)
%==============================================
% Run cutDemo, and kill it when it calls itCut
%
% random walks across scales. set up latent space 
% and test.
%
%  tstLatent.m
%    pickLatent.m
%    setLatent.m
%    coarseFine.m
%      powIter.m
%      gramFixed.m
%      gramSelect.m
%      emKernelFit.m
%
%==============================================
addpath('../SparseDemo');

FALSE = (0 == 1);
TRUE  = ~FALSE;
SANITY  = FALSE;
svdDone = FALSE;

%=========================================
%=== From affinities to markov chains ====
%=========================================
[sA,D,sL,sM] = normalizeAffty(A);
sqrtD    = sparse(diag(D.^0.5));
sqrtDinv = sparse(diag(D .^ -0.5));

%=========================================
%=== Latent space computations        ====
%=========================================
logpow = 2;

% pick kernels for latent space.
latentPar.logpow   = logpow;
latentPar.regular  = 0 ; %regular/irregular (1/0) sampling
latentPar.minRes   = 0.005;
latentPar.boundIdsFlag = 0; % ignore pixels on the boundaries
latentPar.sizeIm   = sizeIm;
[K,selectId] = pickLatent(sM,sqrtD,sqrtDinv,latentPar);
% sort the kernels based on the pixel locations
[selectId,id] = sort(selectId);
% and shuffle the kernels
K = K(:,id);

% set up latent space probabilities (stationary + markov)
u0 = D/sum(D) ;
[Ar,R,st,W] = setLatent(K, u0,1) ;
rbinNhbr = Ar > 0;
rbinNhbr = rbinNhbr - eye(size(rbinNhbr)) ;

%% Debug
% imStats(Ar,Ar')
% figure; plot(sum(R))

% coarse to fine interpolation
[aU,aS] = coarseFine(sL,Ar,K) ;


%===================================================
%===== compare eigenvalues: true vs approximate ====
%===================================================
% svd on the big matrix
[U S V] = svd(full(sL));
%[U S V] = svds(sL,100); 
S = diag(S);

figure(113); clf;
tt  = size(aU,2);
plot(aS(1:tt),'xr-'); hold on;
plot(S(1:tt),'x-'); hold on;

% we can bound the latent space eigenvalues
figure(114); clf;
plot(aS(1:tt),'xr-'); hold on;
plot(S(1:tt).^(2^(logpow-1)),'xb-');
plot(S(1:tt).^(2^(logpow)),'xb-');

% and eigenvectors
for k= 1:size(aU,2)
  figure(111); clf;
  subplot(1,2,1);
  showIm(reshape(U(:,k), sizeIm));
  title(sprintf('Orig Basis: %d (%1.4f)',k,S(k)));
  subplot(1,2,2);
  sgn = sign(U(:,k)' * aU(:,k));
  showIm(reshape(sgn*aU(:,k), sizeIm));
  title(sprintf('Approx Cut: %d (%1.4f)',k,aS(k)));  
  pause;
end


%========================================
%===== kernel intersection          =====
%========================================

% kernel intersection
id1 = 19;
id2 = 15;
pixIds = kernelInter(K(:,id1),K(:,id2),sizeIm,1);

% what the centers of the kernels
%T = zeros(sizeIm);
%T(selectId) = 1;
%figure(114); clf;
%showIm(T);

% kernels with ids drawn as text
%figure(114); clf;
%showIm(im); hold on;
%for k = 1:length(selectId)
%  id = selectId(k);%

%% kernel location
%  idcol = ceil(id/sizeIm(1));
%  idrow = id - (idcol - 1)*sizeIm(1);
%  text(idcol,idrow,sprintf('%d',k),'color','r','fontsize',12,'fontweight','demi');
%end

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

% show neighborhood information for the kernels
for k = 1:length(selectId)

  id = selectId(k);

  figure(115);clf;
  showIm(im); hold on;
  text(kCols,kRows,kText,'color','r','fontsize',12,'fontweight','demi');

  % kernel location
  idcol = ceil(id/sizeIm(1));
  idrow = id - (idcol - 1)*sizeIm(1);
  text(idcol,idrow,kText(k,:),'color','b','fontsize',12,'fontweight','demi');

  % collect neighbors
  nrs = find(rbinNhbr(:,k)) ;
  nrs = selectId(nrs);
  
  % their coordinates
  cols = ceil(nrs/sizeIm(1));
  rows = nrs - (cols - 1)*sizeIm(1);
  
  % mark its neighbors
  nrs = find(rbinNhbr(:,k)) ;
  text(cols,rows,kText(nrs,:),'color','c','fontsize',12,'fontweight','demi');
  
  pause;
end

% find kernel intersections for the links 
% returned in binCut0 (itCut5)
pixIds = [];
debug = 0;
for hostId = 1:size(binCut0,2)
  %pixIds = kernelInter(K(:,id1),K(:,id2),sizeIm,1);
  cutIds = find(binCut0(:,hostId));
  if any(cutIds)
    cutIds = cutIds(cutIds > hostId);
  end
  if any(cutIds)
    for nbrId = 1:length(cutIds)
      if (debug)
	[ hostId cutIds(nbrId) ]
      end
      pixIds = [pixIds kernelInter(K(:,hostId),K(:,cutIds(nbrId)),sizeIm,debug)];
      if (debug)
	pause;
      end
    end
  end
end
% remove duplicated pixel ids
pixIds = unique(pixIds);
pixIds = pixIds(:);
% show pixIds
figure(115);clf;
showIm(im); hold on;
% their coordinates
cols = ceil(pixIds/sizeIm(1));
rows = pixIds - (cols - 1)*sizeIm(1);
%pixIdsText = []
%for t = 1:length(pixIds)
%  pixIdsText = strvcat(pixIdsText,sprintf('%d',pixIds(t)));
%end
%text(cols,rows,pixIdsText,'color','c','fontsize',12,'fontweight','demi');
plot(cols,rows,'cx','linewidth',2);

%=========================================
%=== half-height occupation           ====
%=========================================
[sA,D,sL,sM] = normalizeAffty(A);
sqrtD    = sparse(diag(D.^0.5));
sqrtDinv = sparse(diag(D .^ -0.5));

logpow = 2;
sMp = sM;      
for k = 1:logpow 
  sMp = sMp * sMp;
end
sLp = sqrtDinv * sMp * sqrtD ;  % L = D^-0.5 Markov D^0.5

okIds = 1:size(sLp,2);
boundIdsFlag = 1;
if boundIdsFlag 
  % ids that are not on the boundary
  okIds = boundIds(sizeIm,1);
end

mxM = max(sMp)/2 ;
gtM = sMp > ones(size(sMp,1),1)*mxM;
[mm,okIds] = sort(-mxM);

selectId = [];
pixLoc = [];
done = 0;
k = 1;
while k <= length(okIds) & ~done

  % pick a kernel
  id = okIds(k);
  % pick nhbrs
  nb = find(gtM(:,id));
  % whats common
  co = [];
  co = intersect(pixLoc,nb);
  if ( length(co) < 0.35*length(nb) )
    selectId = [selectId id];
    pixLoc   = [pixLoc nb'];
  end
  
  k = k + 1;
end

pixLoc  = intersect(1:prod(sizeIm),pixLoc);
pixMiss = setdiff(1:prod(sizeIm),pixLoc);

%selectId = sort([selectId pixMiss]);
[mm,tt] = sort(-mxM(pixMiss));
pixMiss = pixMiss(tt);
k = 1;
while k <= length(pixMiss) 

  % pick a kernel
  id = pixMiss(k);
  % pick nhbrs
  nb = find(gtM(:,id));
  % whats common
  co = [];
  co = intersect(pixLoc,nb);
  if ( length(co) < 0.85*length(nb) )
    selectId = [selectId id];
    pixLoc   = [pixLoc nb'];
  end
  
  k = k + 1;
end
pixLoc  = intersect(1:prod(sizeIm),pixLoc);
pixMiss = setdiff(1:prod(sizeIm),pixLoc);

selectId = sort([selectId pixMiss]);

% kernels with ids drawn as text
kCols = ceil(selectId/sizeIm(1));
kCols = kCols(:);
kRows = selectId(:) - (kCols - 1)*sizeIm(1);
kText = [];
for t = 1:length(selectId)
  kText = strvcat(kText,sprintf('%d',selectId(t)));
end
figure(114); clf;
showIm(im); hold on;
text(kCols,kRows,kText,'color','r','fontsize',12,'fontweight','demi');
title(sprintf('# of Kernels: %d',length(selectId)));
