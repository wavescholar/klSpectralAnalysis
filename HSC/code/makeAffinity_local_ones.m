function [Aloc ii jj maxrad] = makeAffinity_local_ones(qaa,minneighbors,maxrad)
% computes affinity function so that each frame is connected to all frames
% within  (qaa-maxrad,qaa+maxrad), and is connected
% to at least a minimum of minneighbors frames (may have to make  more distant
% connections than maxrad)
% edge weight is always 1
%
% n... number of data points. 
% qaa... [n x n_dims] matrix containing location of each data point. This
% is used for finding neighbors
% minneighbors... the minimum number of edges from each point
% maxrad... each point i is connected to all points j within maxrad
% distance in qaa-space. [If there are fewer than minneighbors points in
% this ball, then i is also connected to its closest neighbors, to make at
% least minneighbors neighbors.]
% 
% Local affinities:
%% use sig_i = median(d(rg_i,rg_j)), where j includes all points connected to i.. 
% And
% A(i,j) = 10.^( -d2(si,sj)/(sig_i*sig_j) ),
% See Self-Tuning Spectral Clustering (Zelnik-Manor and Perona)S
% multi version uses e^(sum_k sum_i (rg(i,k) -
% rg(j,k))^2/(sigma_ik*sigma_jk) --> multi-dimensional feature vector rg

nframes = size(qaa,1)

%ndims = size(rg,2); % number of variables 

approxn = 100*minneighbors; % estimate how many neighbors per frame
ii = zeros(approxn*nframes,1);
jj = zeros(approxn*nframes,1);
%ssrg = zeros(approxn*nframes,ndims);

%medrg = zeros(nframes,ndims);

% rgdists = zeros(approxn*nframes,1,'single');

nkeeps = zeros(nframes,1); % record how many neighbors per frame


nij = 0;

for k = 1:nframes
    if mod(k,100) == 1
        disp([num2str(k),': ',num2str(nframes)]);
    end
    
  %  qdist = sum((repmat(qaa(k,:),nframes-k,1)-qaa(k+1:end,:)).^2,2);
    qdist = sum((ones(nframes,1)*qaa(k,:)-qaa).^2,2);
  %  qdist = sum((repmat(qaa(k,:),nframes,1)-qaa).^2,2);
    
    % only look at distances within maxrad
    b = find(qdist <= maxrad);
    
    if size(b,1) < minneighbors+1 % +1 bc includes distance to self
         [~,isort] = sort(qdist);
         b = isort(1:min(minneighbors+1,nframes-k+1));
    end
    
    
   % b = k + b;  % add k, because looking at k+1:end
    
    nb = size(b,1);
    ii(nij + 1:nij + nb) = k*ones(nb,1);
    jj(nij + 1:nij + nb) = b;
  %  ssqdist(nij + 1:nij + nb) = qdist(b);
  
  %  absrgdiff = abs(ones(nb,1)*rg(k,:) - rg(b,:));%  abs(rg(k) - rg(b));
  %  ssrg(nij + 1:nij + nb,:) = absrgdiff;
    
    nij = nij + nb;
    
  %  medrg(k,:) = median(absrgdiff);
    
    nkeeps(k) = nb;
    
end

ii = ii(1:nij);
jj = jj(1:nij);
% ssqdist = ssqdist(1:nij);
%ssrg = ssrg(1:nij,:);

%sigrg = 1.5*medrg;

% if sigrg(end) == 0 % last frame might not have any neighbors other than itself
%     sigrg(end) = 1;
% end
% 
% zerosigrg = find(sigrg == 0);
% if size(zerosigrg,1) > 0
%     sigrg(zerosigrg) = min(sigrg(find(sigrg)));
%     disp('setting zerosigs to min(sigrg(find(sigrg))');
%     size(zerosigrg,1)
% end
% 
% sigrg_ij = sigrg(ii,:).*sigrg(jj,:);


%tt = exp(-(ssrg.^2)./sigrg_ij);
%ttsum = sum(tt,2);

Aloc = sparse(ii,jj,ones(nij,1),nframes,nframes,nij);
%Aloc = sparse(ii,jj,ttsum,nframes,nframes,nij);
% Aloc = sparse(ii,jj,tt,nframes,nframes,nij);



% make sure A is symmetric [sometimes a point j is connected to a point i
% outside of maxrad to meet the minimum number of neighbors, and thus there
% is an edge A(i,j) but not A(j,i). 
u = find(Aloc - Aloc' ~= 0); 
if ~isempty(u) 
    disp('enforcing symmetry');
    
    [ni nj] = ind2sub(size(Aloc),u);
    newvals = zeros(size(ni,1),3);
    nnew = 0;
    for k = 1:size(ni,1)
        if mod(k,1000) == 1
            k
        end
        if Aloc(ni(k),nj(k)) > 0
            if Aloc(nj(k),ni(k)) == 0
                nnew = nnew + 1;
                newvals(nnew,1) = nj(k);
                newvals(nnew,2) = ni(k);
                newvals(nnew,3) = Aloc(ni(k),nj(k));
            end
        end
    end
    newvals = newvals(1:nnew,:);
    nnew
    Asym = sparse(newvals(:,1),newvals(:,2),newvals(:,3),size(Aloc,1),size(Aloc,2),nnew);
    Asym = Asym + Aloc;
    disp(' 0 if symmetric:' );
    max(max(Asym - Asym'))
    
    Aloc = Asym;
end










