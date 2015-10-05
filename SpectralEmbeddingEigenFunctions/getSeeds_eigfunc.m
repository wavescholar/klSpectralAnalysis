function [s m seeds uniqseeds] = getSeeds_eigfunc(Zk,im)

TRUE = 1; FALSE = 0;

dm0 = size(Zk,1);
sizeIm = size(im);

%%%%%%%%%%%% Find seed regions
%% Build latent space, along the lines of Chakra's multi-scale algorithm.
%INSANITY = FALSE;
INSANITY = TRUE;

%% To switch orthogonality testing back on, just change OrthoThresh
%% to some small value (i.e. .2 or.3).

LOOSE_COS_THETA = 0.8;
TIGHT_COS_THETA = 0.985;
MIN_PIXELS = 5;


%%% Plots of Zk
if INSANITY
  
  %% TTD: Fix this to plot just the points in the plane.
  figure(1); clf;
  ide = [1 2 3];
  %ide = [1 3 4];
  % ide = [1 2 6];
  figure(1); hold on;
  plot3(Zk(ide(1),:), Zk(ide(2),:), Zk(ide(3),:),'ob');
  xlabel(int2str(ide(2)));
  ylabel(int2str(ide(3)));
  zlabel(int2str(ide(3)));
  title('Zk_j');
  axis equal;
  hold on;
  
end


%%%%%%%%%%%%%%%%%%%%%%%% Soft clustering of Zk.
%% Sk = Zk/||Zk||
Sk=Zk;                     % same as s in specdiff_mem.m     
Sksum=(sum(Sk.^2,1)).^0.5;	%% Make the vectors unit
Sk=Sk./repmat(Sksum,dm0,1);	%% length

% first pass: tighter control on relative angles
selectId = [];
pixLoc   = ones(prod(sizeIm),1)>0;
nb = zeros(prod(sizeIm),1);
okIds = find(pixLoc);
while length(okIds) > 0
  
  % pick a kernel
  id=floor(rand*(length(okIds)-1))+1;
  x = okIds(id);
  selectId = [selectId x];
  Skx = Sk(:,x);
  
  % delete kernel and nhbrs from pixLoc.
  hit_matrix = Skx'*Sk;  %% Calculate the cos(theta) with all other
  
  % doing the comparison where necessary 
  pixLoc(hit_matrix > LOOSE_COS_THETA) = FALSE;
  
  okIds = find(pixLoc);
end

nCenter = length(selectId);

%% Assign pixels to centers.
centers = Sk(:,selectId); % same as m in zeg.m (except possibly made differently)
hit_matrix = centers'*Sk; 
if size(hit_matrix,1) > 1
  [mx mxHit] = max(hit_matrix);
else
  mxHit = ones(1, size(hit_matrix,2));
  mx = hit_matrix;
end
centerLabel = mxHit;
centerLabel(mx <= TIGHT_COS_THETA) = 0;

if INSANITY
  
  figure(1); clf;
  ide = [1 2 3];
  %ide = [1 3 4];
  % ide = [1 2 6];
  figure(1); hold on;
%  plot3(Zk(ide(1),:), Zk(ide(2),:), Zk(ide(3),:),'ok');
  plot(1000*Zk(ide(2),:), 1000*Zk(ide(3),:),'o','Color',[.35 .35 .35]);
%  xlabel(int2str(ide(2)));
%  ylabel(int2str(ide(3)));
  xlabel('z_2 x 1000');
  ylabel('z_3 x 1000');
%  zlabel(int2str(ide(1)));
%  title('Zk_j');
  axis equal;axis([-25 20 -25 20]);
  hold on;
  
  %% Display the center initializations in Z.
  RandCols=(.8*rand(nCenter,3))+.2;
  for j=1:nCenter
    idx = (centerLabel == j);
    figure(1); hold on;
%    plot3(Zk(ide(1),idx), Zk(ide(2),idx), Zk(ide(3),idx),'*', 'Color', ...
%          RandCols(j,:));
    plot(1000*Zk(ide(2),idx), 1000*Zk(ide(3),idx),'*', 'Color', ...
         RandCols(j,:));
    ii=find(idx>0);
    h=text(1000*Zk(ide(2),ii(1)), 1000*Zk(ide(3),ii(1)), int2str(j), 'FontSize', 16, 'Color', [0 0 0]);
  end  
  
  %% Display the center initializations in an image.
  centerIm=zeros(sizeIm(1)*sizeIm(2),3);
  for j=1:nCenter
    idx = (centerLabel == j);
    centerIm(idx,:)=repmat(RandCols(j,:),sum(idx),1);
  end;

  figure(700);
  clf;
  image(reshape(centerIm,sizeIm(1),sizeIm(2),3).*repmat((im/512)+.5,[1 1 3]));
  axis equal; axis off;
  title('Initial cluster centroids');

end

%% Second pass: Adjust these center locations according to weighted
%% means.
for k = 1:3
  centerLabel = mxHit;
  centerLabel(mx <= LOOSE_COS_THETA) = 0;

  sumWght = zeros(nCenter,1);
  newCenter = zeros(dm0, nCenter);
  wght = hit_matrix - LOOSE_COS_THETA;
  wght = max(wght,0);
  for j=1:nCenter
    idx = centerLabel(:)==j;
    wght(j,~idx) = 0;
    sumWght(j) = sum(wght(j,idx));
    if sumWght(j) > 0
      newCenter(:,j) = sum(repmat(wght(j,idx), dm0,1) .* Zk(:,idx), 2)/ ...
          sumWght(j);
      newCenter(:,j) = newCenter(:,j)/sqrt(sum(newCenter(:,j).^2));
    end
  end
  hit_matrix = newCenter'*Sk; 
  if size(hit_matrix,1) > 1
    [mx mxHit] = max(hit_matrix);
  else
    mxHit = ones(1, size(hit_matrix,2));
    mx = hit_matrix;
  end

  if INSANITY
    %% Display the center initializations in an image.
    centerLabel = mxHit;
    centerLabel(mx <= TIGHT_COS_THETA) = 0;
    centerIm=zeros(sizeIm(1)*sizeIm(2),3);
    for j=1:nCenter
      idx = (centerLabel == j);
      centerIm(idx,:)=repmat(RandCols(j,:),sum(idx),1);
    end;

    figure(700+k);
    clf;
    image(reshape(centerIm,sizeIm(1),sizeIm(2),3).*repmat((im/512)+.5,[1 1 3]));
    axis equal; axis off;
    title('Initial cluster centroids');
  end
end

%% Third pass: Remove any centers with zero weight. Add centers to cover 
%% uncovered pixels.
centerLabel = mxHit;
centerLabel(mx <= TIGHT_COS_THETA) = 0;
for j=1:size(newCenter,2)
  sumWght(j) = sum(centerLabel==j);
end
centers = newCenter(:, sumWght>0);
nCenter = size(centers,2);

hit_matrix = centers'*Sk; 
if size(hit_matrix,1) > 1
  [mx mxHit] = max(hit_matrix);
else
  mxHit = ones(1, size(hit_matrix,2));
  mx = hit_matrix;
end
uncovered = find(mx < LOOSE_COS_THETA);

while length(uncovered)>0
  
  % pick a kernel
  id=floor(rand*(length(uncovered)-1))+1;
  x = uncovered(id);
  
  %% compute inner-products
  centers = [centers Sk(:,x)];
  hit_matrix = [hit_matrix; Sk(:,x)' * Sk];
  
  %% compute max assignments
  if size(hit_matrix,1) > 1
    [mx mxHit] = max(hit_matrix);
  else
    mxHit = ones(1, size(hit_matrix,2));
    mx = hit_matrix;
  end
  uncovered = find(mx < LOOSE_COS_THETA);
  
end

%% Set variables used by runCut...
centerLabel = mxHit;
centerLabel(mx <= TIGHT_COS_THETA) = 0;
nSeedReg = size(centers,2);
seedReg = zeros(prod(sizeIm),nSeedReg)>0;
for k=1:nSeedReg
  seedReg(centerLabel == k,k) = TRUE;
end

%% Throw away any seeds with less than MIN_PIXELS pixels.
if MIN_PIXELS > 0
  %% Small regions
  MIN_PIXELS = 5;
  seedReg_old=seedReg;
  centers_old=centers;

  sInd=1;
  maxRegSize=0;

  clear seedReg0;
  clear centers0;

  for i=1:1:size(seedReg,2)
    if (length(find(seedReg(:,i)>0))>MIN_PIXELS)
      seedReg0(:,sInd)=seedReg(:,i);
      centers0(:,sInd)=centers(:,i);
      sInd=sInd+1;
    end;
  end;


  if size(seedReg0,2) ~= size(seedReg,2)
    %% Recompute hit matrix for new centers
    hit_matrix = centers0'*Sk; 
    if size(hit_matrix,1) > 1
      [mx mxHit] = max(hit_matrix);
    else
      mxHit = ones(1, size(hit_matrix,2));
      mx = hit_matrix;
    end
    
    %% Set variables used by runCut...
    centerLabel = mxHit;
    centerLabel(mx <= TIGHT_COS_THETA) = 0;
    nSeedReg = size(centers,2);
    seedReg = zeros(prod(sizeIm),nSeedReg)>0;
    for k=1:nSeedReg
      seedReg(centerLabel == k,k) = TRUE;
    end
    seedReg=seedReg0;
    centers=centers0;
  end
end

fprintf(2,'Number of seeds: %d\n', nSeedReg);

if ~exist('RandCols','var') | size(RandCols,1) ~= size(seedReg,2)
  RandCols=(.8*rand(size(seedReg,2),3))+.2;
end


% vmb:
rands = rand(size(seedReg,2),3);
figure; imagesc(im); hold on;
Seeds = cell(size(seedReg,2),1);
for j = 1:size(seedReg,2)
    Seeds{j} = find(seedReg(:,j));
    [a b] = ind2sub(size(im),Seeds{j});
    plot(b,a,'ms','MarkerSize',10,'MarkerFaceColor',rands(j,:));
end
title('Final cluster seeds','FontSize',20);

%% adjust Seeds so that all cells have same number of entries
M = size(seedReg,2);
nseeds = zeros(M,1);
for k = 1:M
    nseeds(k) = size(Seeds{k},1);
end
Nseeds = max(nseeds);
seeds = cell(size(Seeds));
for k = 1:M
    seeds{k} = zeros(Nseeds,1);
    seeds{k}(1:nseeds(k)) = Seeds{k};
end
seeds = cell2mat(seeds');
uniqseeds = seeds;

m = centers; 
s = Sk;



%% See what's left!
centerIm=zeros(sizeIm(1)*sizeIm(2),3);
for j=1:size(seedReg,2);
  idx = seedReg(:,j);
  centerIm(idx,:)=repmat(RandCols(j,:),sum(idx),1);
end;
 figure(701);
 clf;
 %image(reshape(centerIm,sizeIm(1),sizeIm(2),3).*repmat((im/512)+.5,[1 1 3]));
 centerIm(centerIm==0)=.35;
 image(reshape(centerIm,sizeIm(1),sizeIm(2),3));
 [x y] = meshgrid(1:sizeIm(2), 1:sizeIm(1));
 for j=1:size(seedReg,2);
   mxy = [sum(x(:) .*seedReg(:,j)), ...
          sum(y(:) .*seedReg(:,j))]/sum(seedReg(:,j));
   if sum(RandCols(j,:))/3> 0.7
     h=text(mxy(1), mxy(2), int2str(j), 'FontSize', 16, 'Color', [.01 .01 .01]);
   else
     h=text(mxy(1), mxy(2), int2str(j), 'FontSize', 16, 'Color', [1 1 1]);
   end      
 end;
 axis equal; axis off;
 %title('Cluster centroids');

  
