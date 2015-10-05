function [Aloc ssfeats medfeats ii jj maxrad] = makeAffinity_local_cerno(feats,minneighbors,maxrad)
% computes affinity function so that each frame is connected to all frames
% within  (data-maxrad,data+maxrad), and is connected
% to at least a minimum of minneighbors frames (may have to make  more distant
% connections than maxrad)
%
% n... number of data points.
% data... [n x n_dims] matrix containing location of each data point. This
% is used for finding neighbors [ REPLACED BY FEATS for cerno. ]
% feats... a vector/matrix with values at each point [could be intensity,
% color,...]. n * n_features. This is used in the affinity function
% minneighbors... the minimum number of edges from each point
% maxrad... each point i is connected to all points j within maxrad
% distance in data-space. [If there are fewer than minneighbors points in
% this ball, then i is also connected to its closest neighbors, to make at
% least minneighbors neighbors.]
%
% Local affinities:
%% use sig_i = median(d(feats_i,feats_j)), where j includes all points connected to i..
% And
% A(i,j) = 10.^( -d2(si,sj)/(sig_i*sig_j) ),
% See Self-Tuning Spectral Clustering (Zelnik-Manor and Perona)S
% multi version uses e^(sum_k sum_i (feats(i,k) -
% feats(j,k))^2/(sigma_ik*sigma_jk) --> multi-dimensional feature vector feats

debug = 0; 

nframes = size(feats,1)

ndims = size(feats,2); % number of variables
nfeats = size(feats,2);

approxn = 10*minneighbors; % estimate how many neighbors per frame
ii = zeros(approxn*nframes,1);
jj = zeros(approxn*nframes,1);
%ssabsfeats = sparse(approxn*nframes,ndims);
ssfeats = zeros(approxn*nframes,ndims,'single');

medfeats = zeros(nframes,ndims);
mu_hat = zeros(nframes,ndims);
b_hat = zeros(nframes,ndims);

% featsdists = zeros(approxn*nframes,1,'single');

nkeeps = zeros(nframes,1); % record how many neighbors per frame


nij = 0;

for k = 1:nframes
        if mod(k,100) == 1
                disp([num2str(k),': ',num2str(nframes)]);
        end
        
        %  qdist = sum((repmat(data(k,:),nframes-k,1)-data(k+1:end,:)).^2,2);
        qdist = sum((ones(nframes,1)*feats(k,:)-feats).^2,2);
        %  qdist = sum((repmat(data(k,:),nframes,1)-data).^2,2);
        
        % only look at distances within maxrad
        b = find(qdist <= maxrad);
        
        if size(b,1) < minneighbors+1 % +1 bc includes distance to self
                [uuu,isort] = sort(qdist); clear uuu;
                b = isort(1:min(minneighbors+1,nframes-k+1));
        end
        
        
        % b = k + b;  % add k, because looking at k+1:end
        
        nb = size(b,1);
        ii(nij + 1:nij + nb) = k*ones(nb,1);
        jj(nij + 1:nij + nb) = b;
        %  ssqdist(nij + 1:nij + nb) = qdist(b);
        
       % absfeatsdiff = abs(ones(nb,1)*feats(k,:) - feats(b,:));%  abs(feats(k) - feats(b));
        featsdiff = ones(nb,1)*feats(k,:) - feats(b,:);
    %    ssabsfeats(nij + 1:nij + nb,:) = abs(featsdiff);
        
        % prepare for laplacian:
        ssfeats(nij + 1:nij + nb,:) = featsdiff;
        mu_hat(k,:) = median(featsdiff);
        b_hat(k,:) = sum(abs(featsdiff - ones(nb,1)*mu_hat(k,:)))/nb;
        
        nij = nij + nb;
        
        medfeats(k,:) = median(abs(featsdiff));
        
        nkeeps(k) = nb;
        
end

ii = ii(1:nij);
jj = jj(1:nij);
ssfeats = ssfeats(1:nij,:);

if debug
        nbins = 50;
        for f = 1:nfeats
                
                mindists = ssfeats(:,f);
                h = figure;
                [n s] = hist(mindists(:),nbins);
                n = n/sum(n);
                semilogy(s,n,'r','LineWidth',2);
                
                formatfig;
                % title('Distribution of (x - y) for each normalized (mu 0 std 1) feature. Fit to gaus (green), lapl (blue)');
                title(num2str(f));
                
                hold on;
                
                % Laplacian:
                % generate paramaters mu_hat and b_hat:
                muhat = median(mindists);
                bhat = sum(abs(mindists - muhat))/length(mindists);
                % approximate distribution with Laplacian:
                lap = @(x) (1/(2*bhat))*exp(-(abs(x - muhat)/bhat));
                %   plot(mindists,lap(mindists),'b.','LineWidth',2);
                la = lap(mindists);
                %    la = la/max(la)*max(n);
                semilogy(mindists,la,'b.','MarkerSize',2);
                
                % Gaussian:
                % generate paramaters mu_hat and b_hat:
                muhat = mean(mindists);
                stdhat = std(mindists);
                % approximate distribution with gaussian:
                gauss = @(x) (1/(stdhat*sqrt(2*pi)))*exp(-((x-muhat).^2/(2*stdhat^2)));
                %     plot(mindists,gauss(mindists),'y.','LineWidth',2);
                g = gauss(mindists);
                %       g = g/max(g)*max(n);
                semilogy(mindists,g,'g.','MarkerSize',2);
                
                %         % Exponential:
                %         % generate paramaters mu_hat and b_hat:
                %         mu_hat = mean(mindists);
                %         lambda_hat = 1/mu_hat;
                %         % approximate distribution with exponential:
                %        % need if, then, so define an actual function
                %         hold on;  plot(mindists,expo(mindists,lambda_hat),'g.','LineWidth',2);
                %
                %         % Chi-squared
                %          k = 2;
                %         hold on;  plot(mindists,chisquared(mindists,k),'m.','LineWidth',2);
                % saveas(h,[figpath,'normed_(x-y)_feature_',num2str(f),'_distributions_compare_approx_nrand',num2str(nrand),'_minneighbors',num2str(minneighbors),'.jpg'])
                pause;
                close all;
        end
end

% ssqdist = ssqdist(1:nij);
%ssabsfeats = abs(ssfeats(1:nij,:));

if 0 % gaussian distribution:
        sigfeats = 1.5*medfeats;
        
        for k = 1:nfeats
                if sigfeats(end,k) == 0
                        sigfeats(end,k) = 1;
                end
        end
        
        zerosigfeats = find(sigfeats == 0);
        if size(zerosigfeats,1) > 0
                sigfeats(zerosigfeats) = min(sigfeats(find(sigfeats)));
                disp('setting zerosigs to min(sigfeats(find(sigfeats))');
                size(zerosigfeats,1)
        end
        
        sigfeats_ij = sigfeats(ii,:).*sigfeats(jj,:);
        
        tt = exp(-(abs(ssfeats(1:nij,:)).^2)./sigfeats_ij);
end

for k = 1:nfeats
        if b_hat(end,k) == 0
                b_hat(end,k) = 1;
        end
end

zerobhatfeats = find(b_hat == 0);
if size(zerobhatfeats,1) > 0
        b_hat(zerobhatfeats) = min(b_hat(find(b_hat)));
        disp('setting zero b_hats to min(b_hats)) - already corrected for final tile');
        size(zerobhatfeats,1)
end


% 
bhats_ij = sqrt(b_hat(ii,:).*b_hat(jj,:));

tt = (1./(2*bhats_ij)).*exp(-abs(ssfeats)./bhats_ij   );
ttsum = sum(tt,2);
ttsum = double(ttsum/nfeats);

%
% tt = -(ssabsfeats.^2)./sigfeats_ij;
% ttsum = sum(tt,2);
% ttsum = exp(ttsum);

Aloc = sparse(ii,jj,ttsum,nframes,nframes,nij);
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

