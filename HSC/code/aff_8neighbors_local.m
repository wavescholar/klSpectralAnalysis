function [A,sigma,dists,ii,jj] = aff_8neighbors_local(im,findxx,zthresh,windows)
% gets sparse affinity matrix for image, looking at distances to 8-neighbors of
% each pixel, and only looking at pixels in findxx
% a(i,j) = exp(-|xi - xj|^2/(sig_i*sig_j)), for i and j within 8-neighborhood of i, 
% and sig_i = std(|xi - xj|) for for i and j within 8-neighborhood of i (analogue
% for j)
% computes SIGMA = median(abs(dists))/2

if nargin < 3
    zthresh = 0; % only look at neighbors (in 8-window) > zthresh
end
if nargin < 4 % windows = size of window around each point (window = 1: 8 neighbors, window = 2: 24 neighbors,...
    windows = 1;
end

figs = 0;

w = (1+2*windows)^2 - 1; % maximum number of pixels in window
n = size(findxx,1); % number of pixels in findxx

maxn = w*n; % maximum number of distances (to pre-allocate matrix)

% add zero-border to image (so that it is easier to deal with pixels on
% border of image)
[nrows ncols] = size(im); % number of rows and columns in image
[fx fy] = ind2sub([nrows ncols],findxx);
Im = zeros(nrows + 2*windows, ncols + 2*windows);
Im(windows+1:end-windows,windows+1:end-windows) = im;
newfx = fx + windows; newfy = fy + windows;
[newnrows newncols] = size(Im); 
newfindxx = sub2ind([newnrows newncols],newfx,newfy);
% check:
min(Im(newfindxx) == im(findxx))

im = Im;
findxx = newfindxx;
nrows = newnrows; ncols = newncols;

% make a look-up table from Im to findxx (wrt to new zero-bordered image)
% - used for creating sparse matrix
lookup = zeros(prod(size(im)),1);
for k = 1:n
    lookup(findxx(k)) = k;
end

% now get affinities:
ii = zeros(maxn,1);
jj = zeros(maxn,1);
dists = zeros(maxn,1);
sigma = zeros(n,1); % compute sigma locally for each pixel
tot = 0;

for k = 1:n
    
    if mod(k,1000) == 0
        disp([num2str(k),': ',num2str(n)]);
    end
    t = findxx(k);
    % get neighborhood of t:
    
    neighbors = zeros(1,w);
    nt = 0;
    for nc = -windows : windows
        ncc = t + nc*nrows;
        for nr = -windows : windows
            nt = nt + 1;
            neighbors(nt) = ncc + nr;
        end
    end
    
     % find non-zero neighbors
    imneigh = im(neighbors);
    findneigh = find(imneigh > zthresh);
    neighbors = neighbors(findneigh);
    imneigh = imneigh(findneigh);
    
    if figs
        % check:
        figurev(110); clf;
        imagesc(im); colormap(gray);
        hold on;
        [a b] = ind2sub(size(im),t);
        plot(b,a,'m*','MarkerSize',10);
        [a b] = ind2sub(size(im),neighbors);
        plot(b,a,'gx','MarkerSize',5);
        pause(.005);
    end
    
    % record new entries in A:
    nneigh = size(imneigh,2);
    ii(tot + 1:tot + nneigh) = repmat(k,nneigh,1);
    jj(tot + 1:tot + nneigh) = lookup(neighbors); % index of neighbors in findxx
    distsk = im(t) - im(neighbors);
    dists(tot + 1:tot + nneigh) = distsk;
    sigma(k) = std(abs(distsk));
    if sigma(k) == 0
        sigma(k) = 1;
    end

    tot = tot + nneigh;
end

ii = ii(1:tot);
jj = jj(1:tot);
dists = dists(1:tot);

% find sigma_i*sigma_j for each pair
sigij = sigma(ii).*sigma(jj);

  %  SIGMA = 1.5*median(abs(dists));

%disp(['SIGMA = ',num2str(SIGMA)]);

% ss = exp(-0.5*dists.^2/SIGMA^2);
ss = exp(-dists.^2./sigij);

A = sparse(ii,jj,ss,n,n,tot);


