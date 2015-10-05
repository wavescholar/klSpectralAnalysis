%% code for hierarchical network segmentation:

addpath('svdms');
addpath('svdms/Latent');
addpath('code');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex 1: segment data points in 2D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load testdata1;
figure; scatter(data(:,1),data(:,2),[],vals,'filled');

% Create a connectivity matrix (A) for this data:

nrand = 1000;
minneighbors = 5; 

nframes = size(data,1)
rr = floor(nframes*rand(nrand,1)) + 1;
mindists = zeros(nrand,1);
for k = 1:nrand
    kk = rr(k);
    %    qdist = sum((repmat(qaa(kk,:),nframes,1)-qaa(1:nframes,:)).^2,2);
    qdist = sum((ones(nframes,1)*data(kk,:)-data).^2,2);
    qdists = sort(qdist);
    mindists(k) = qdists(minneighbors);
end
maxrad = median(mindists)

[A ssrg medrg ii jj maxrad] = makeAffinity_local(data,vals,minneighbors,maxrad);

% take a look at the affinity matrix: 
figure; spy(A,'k'); title('A','FontSize',20); set(gca,'FontSize',12); axis square;


%% HIERACHICAL NETWORK SEGMENTATION!!!
subsp = 2;
nrounds = 50;
hh = svdms_msm(A,subsp,nrounds);
% Interesting: in svdms_msm.m, comment out return on line 301 to compute
% eigenvalues and eigenvectors.

if size(hh{end}.A,1) == 1
    hh = hh(1:end-1);
end

% visualize affinity matrix at each hierarchy level:
nhh = size(hh,2);
if floor(sqrt(nhh)) == sqrt(nhh)
    ns = sqrt(nhh)
else
    ns = floor(sqrt(nhh)) + 1
end
figure;
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k');
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end


Colors = {'y';'m';'c';'r';'g';'b';'k'};

% visualize segmentation at each hierarchy level
figure;
for lev = 2:nhh
    subplot(ns,ns,lev);
    if max(hh{lev}.cl) > 7
    scatter(data(:,1),data(:,2),[],hh{lev}.cl,'filled');
    else
        hold on;
        for k = 1:max(hh{lev}.cl)
            u = find(hh{lev}.cl == k);
            scatter(data(u,1),data(u,2),[],Colors{k},'filled');
        end
    end
    title(['level ',num2str(lev),', ',num2str(size(hh{lev}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12); 
end

% visualize segmentation in affinity map... 
figure;
for lev = 2:nhh
    hh{lev}.colA = zeros(size(A));
    hh{lev}.colA(find(A)) = .5;
    subplot(ns,ns,lev); 
    if max(hh{lev}.cl) > 7
        for k = 1:max(hh{lev}.cl)
            u = find(hh{lev}.cl == k);
            hh{lev}.colA(u,u) = k;
        end
        imagesc(hh{lev}.colA);
    else
        
        hold on;
        for k = 1:max(hh{lev}.cl);
            u = find(hh{lev}.cl == k);
            Au = zeros(size(A));
            Au(u,u) = 1;
            spy(Au,Colors{k});
        end
    end
    title(['level ',num2str(lev),', ',num2str(size(hh{lev}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex 2: segment data points in 3D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load testdata2;
figure; scatter3(data(:,1),data(:,2),data(:,3),[],vals,'filled');

% Create a connectivity matrix (A) for this data:

nrand = 500;

minneighbors = 5; 

nframes = size(data,1)
rr = floor(nframes*rand(nrand,1)) + 1;
mindists = zeros(nrand,1);
for k = 1:nrand
    kk = rr(k);
    %    qdist = sum((repmat(qaa(kk,:),nframes,1)-qaa(1:nframes,:)).^2,2);
    qdist = sum((ones(nframes,1)*data(kk,:)-data).^2,2);
    qdists = sort(qdist);
    mindists(k) = qdists(minneighbors);
end
maxrad = median(mindists)

[A ssrg medrg ii jj maxrad] = makeAffinity_local(data,vals,minneighbors,maxrad);

% take a look at the affinity matrix: 
figure; spy(A,'k'); title('A','FontSize',20); set(gca,'FontSize',12); axis square;


%% HIERACHICAL NETWORK SEGMENTATION!!!
subsp = 2;
nrounds = 50;
hh = svdms_msm(A,subsp,nrounds);

if size(hh{end}.A,1) == 1
    hh = hh(1:end-1);
end

% visualize affinity matrix at each hierarchy level:
nhh = size(hh,2);
nhh = size(hh,2);
if floor(sqrt(nhh)) == sqrt(nhh)
    ns = sqrt(nhh)
else
    ns = floor(sqrt(nhh)) + 1
end
figure;
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k');
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end


Colors = {'y';'m';'c';'r';'g';'b';'k'};

% visualize segmentation at each hierarchy level
figure;
for lev = 2:nhh
    subplot(ns,ns,lev);
    if max(hh{lev}.cl) > 7
    scatter3(data(:,1),data(:,2),data(:,3),[],hh{lev}.cl);
    else
        hold on;
        for k = 1:max(hh{lev}.cl)
            u = find(hh{lev}.cl == k);
            scatter3(data(u,1),data(u,2),data(u,3),[],Colors{k});
        end
    end
    title(['level ',num2str(lev),', ',num2str(size(hh{lev}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end

% visualize segmentation in affinity map... 
figure;
for lev = 2:nhh
    hh{lev}.colA = zeros(size(A));
    hh{lev}.colA(find(A)) = .5;
    subplot(ns,ns,lev); 
    if max(hh{lev}.cl) > 7
        for k = 1:max(hh{lev}.cl)
            u = find(hh{lev}.cl == k);
            hh{lev}.colA(u,u) = k;
        end
        imagesc(hh{lev}.colA);
    else
        
        hold on;
        for k = 1:max(hh{lev}.cl);
            u = find(hh{lev}.cl == k);
            Au = zeros(size(A));
            Au(u,u) = 1;
            spy(Au,Colors{k});
        end
    end
    title(['level ',num2str(lev),', ',num2str(size(hh{lev}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex 3: segment an image:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load testim;
im = double(im); 
figure; imagesc(im); colormap(gray); axis off; title('test image','FontSize',20);


% The example here uses an image. The only real differences between using images
% and data points is how the affinity matrix is defined (the affinity 
% matrix used here for images is a special case of the matrix used for data points)
% and how the images
% are plotted (use imagesc for images, scatter for data points). 
% See example 1-2 below using data points in space!

% Create a connectivity matrix (A) for this image:
zthresh = 0;
findxx = find(im > zthresh);
windows = 1; % standard 8 neighborhood
[A,sigma,dists,ii,jj] = aff_8neighbors_local(im,findxx,zthresh,windows);

% similarly, could create a connectivity matrix bw points in a network
% [A ssrg medrg ii jj maxrad] = makeAffinity_local(qaa,rg,minneighbors,maxrad);


% take a look at the affinity matrix: 
figure; spy(A,'k'); title('A','FontSize',20); set(gca,'FontSize',12);


%% HIERACHICAL NETWORK SEGMENTATION!!!
subsp = 2;
nrounds = 50;
hh = svdms_msm(A,subsp,nrounds);

if size(hh{end}.A,1) == 1
    hh = hh(1:end-1);
end

% visualize affinity matrix at each hierarchy level:
nhh = size(hh,2);
nhh = size(hh,2);
if floor(sqrt(nhh)) == sqrt(nhh)
    ns = sqrt(nhh)
else
    ns = floor(sqrt(nhh)) + 1
end
figure;
for k = 1:nhh
    subplot(ns,ns,k);
    spy(hh{k}.A,'k');
    title(['level ',num2str(k),', ',num2str(size(hh{k}.A,1)),' clusters'],'FontSize',16);
    set(gca,'FontSize',12);
end

% visualize segmentation at each hierarchy level
figure;
for lev = 2:nhh
    subplot(ns,ns,lev);
    clim = zeros(size(im));
    clim(findxx) = hh{lev}.cl;
    imagesc(clim); colormap(gray);
    title(['level ',num2str(lev),', ',num2str(size(hh{lev}.A,1)),' clusters'],'FontSize',16);
    axis image;
    set(gca,'FontSize',12);
end





%% something else:
% pertubation of affinity matrix to group clusters based on sensitivity of
% edges:

lev = 2;
sizeim = size(im);
%if lev > 1
    [hh{lev}.Acut hh{lev}.bottlenecks] = perturb_easy_im(hh{lev},findxx,sizeim);
% else % run k-means
%     [hh{lev}.Acut hh{lev}.bottlenecks] = perturb_easy_im(hh{lev},findxx,sizeim,im);
% end

tau = 0;
[Acc idividefurther.indices] = cca_wrap(hh{lev}.Acut,tau);
assn = zeros(size(hh{lev}.Acut,1),1);
nccacl = size(idividefurther.indices,1);
for k = 1:nccacl
    assn(idividefurther.indices{k}) = k;
end

clsassn = assn(hh{lev}.cl);
figure; 
clim = zeros(size(im));
clim(findxx) = clsassn;
imagesc(clim);

% remove clusters that have low intensities
ncl = max(clim(:));
meanint = zeros(ncl,1);
for k = 1:ncl
    u = find(clim == k);
    meanint(k) = mean(im(u));
end
climint = meanint(clim);
figure; imagesc(climint); colorbar;
figure; plot(meanint);

% for now: handpick cutoff - very sketchy
cutoff =  climint(10,8)%climint(12,21)%climint(1,1); 
brightcls = find(meanint >= cutoff);
darkcls = find(meanint < cutoff);
idarkcls = [];
for k = 1:size(darkcls,1)
    idarkcls = [idarkcls; find(clim == darkcls(k))];
end

climint_bright = climint;
climint_bright(idarkcls) = 0;
clim_bright = clim;
clim_bright(idarkcls) = 0;

figure; imagesc(climint_bright); 
figure; imagesc(clim_bright); 

ncl = size(unique(clim_bright),1) - 1; % leave out zero

% remove tiny clusters
cls = cell(ncl,1);
sizecls = zeros(ncl,1);
for k = 1:ncl
    kb = brightcls(k);
    cls{k} = find(clim_bright == kb);
    sizecls(k) = size(cls{k},1);
end

tinycutoff = 5;
keepcls = brightcls(find(sizecls > tinycutoff));

% renumber clusters
clim2 = zeros(size(clim_bright));
ncl = size(keepcls,1);
cls = cell(ncl,1);
clsizes = zeros(ncl,1);
for k = 1:ncl
    kb = keepcls(k);
    u = find(clim_bright == kb);
    clim2(u) = k;
    clsizes(k) = size(u,1);
    cls{k} = u;
end
clim = clim2; clear clim2;

% now change blobs to outlines!!!
% and make results images
bwclim = zeros(size(clim));
bwclim(find(clim)) = 1;
NB = bwboundaries(bwclim,8);

clborders = im;
figure; imagesc(clborders); colormap(gray);
hold on;
for k=1:length(NB)
    border = NB{k};
    plot(border(:,2),border(:,1),'rs','MarkerFaceColor','r','MarkerSize',10);
end





