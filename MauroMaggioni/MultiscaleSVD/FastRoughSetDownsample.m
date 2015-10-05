function vNet = FastRoughSetDownsample( cX, cOpts )

%
% function vNet = FastRoughSetDownsample( cX, cOpts )
%
% IN:
%   [cX]        : D by N matrix of N points in D dimensions.
%   cOpts       : stricture with options
%                   Delta               : find a maximal cOpts.Delta net. If cOpts.kNN is not empty, this is disregarded
%                   kNN                 : find a maximal cOpts.kNN net. 
%                   [Dim]               : dimensionality of the data. Needed if cX is not provided. Default: dimension of cX.
%                   [RestrictSetIdxs]   : restricts the net to points with these indices. Default: 1:M (all points allowed)
%                   [Sigma]             : exponential weights for each chosen point: e^(-(dist/cOpts.Sigma)^2). Default: Delta
%                   [DownSample]        : downsamples the set around the net. Default: true.
%                   [DoSVD]             : compute local svd. Default: false.
%                   [SVDDim]            : how many dimensions to compute in the SVD. If in metric space, it is the dimension of the embedding space
%                                         for multidimensional scaling. Default: Dim.
%                   [DoD]               : compute local weight D (see outputs). Default: false.
%                   [SVDDelta]          : compute SVD on points in a SVDDelta-ball around each point.
%                                         SVDDelta should be >=Delta. Default: cOpts.Delta.
%                   [NormalizedSVD]: 0,1. Normalizes the local singular values by the number of local points
%                                         used to compute it. Default: 1.
%                   [RandomizedSVD]     : uses Rokhlin-Tygert randomized pca algorithm. Default: false.
%                   [KeepV]             : whether to keep the principal components or not (memory intensive). Default: 1.
%                   [KeepCenter]        : whether to keep the center of each cell in the net. Default: 0.
%                   [AvoidSmallClusters]: re-assigns points in clusters of size smaller than AvoidSmallClusters*(median cluster size) to the nearest clusters.
%                                         Does not update the SVD data nor D for the enlarged clusters. Default: 0 (no reassignment).
%                   [Randomized]        : randomizes indices of the points. Default : 1.
%                   [DistInfo]          : if not empty, it will be used for distance computation. It is assumed that it contains at least
%                                        the neighbors up to the requested distance Radius (and PCARadius if provided), or kNN (if provided).
%                                        It assumes the neighbors are listed in increasing distance. It's a structure containing the following fields:
%                           count           : N vector, the k-th entry is the number of neighbors of the n-th point
%                           idxs            : N cell, the k-th cell contains the indices of the nearest neighbors of the k-th point
%                           dists           : N cell, the k-th cell contains the distances of the nearest neighbors of the k-th point (corresponding to idxs).
%                                             Should be sorted in incrasing magnitude.
%                   [RandomizedSVDThres]: if the number of points in a cOpts.Delta neighborhood exceeds this, use randomized svd computation. Default: Inf.
%                   [NNInfo]            : extra structure for nearest neighbor extra information. May contain:
%                       Atria           : If the TSTOOL nn_prepare is used, this is the atria structure of cX. It will be used for faster computation.
%                   [Verbose]           : 0,1. Default: 1.
%
%
% OUT:
%   vNet.Idxs       : cX(:,vNet.Idxs) is the maximal cOpts.Delta net
%   vNet.InvIdxs    : vNet.InvIdxs is an M column vector, and cX(:,vNet.Idxs(vNet.InvIdxs(j))) is the cluster center point to which
%                       the point cX(:,j) is assigned
%   vNet.S          : |vNet.Idxs| by D matrix of singular values of local covariance matrix at each net point. Empty if ~cOpts.DoSVD
%   vNet.V          : |vNet.Idxs| by D by D matrix of singular vectors of local svNet.D at each point. Empty if ~cOpts.DoSVD
%   vNet.D          : vNet.D(k) is the mass to assign to the k-th selected point cX(:,vNet.Idxs(k)). Does not include mass of point itself.
%   vNet.Center     : |vNet.Idxs| by D matrix of the center of each cell in the D dim space. Not returned if  ~KeepCenter
%   vNet.Count      : number of neighbors around each point used for the local svd computation
%   vNet.NNInfo     : extra info for nearest neighbors, as returned by nrsearch
%   vNet.Delta      : delta used for the net construction. If a nearest neighbor net was constructed, this is an N vector, otherwise a scalar.
%
% USES:
%   nn_prepare, nn_search, range_search         : from OpenTSTool package, if available
%   ANNsearch,                                  : ANN package by D. Mount, if available
%   nrsearch,
%   setdiff_sorted                              : from LightSpeed package by T. Minka
%   pca                                         : by M. Tygert
%
%
% EXAMPLE:
%
%   See also FastRoughSetDownsample_Test for examples, and tests of computational complexity.
%
%
% Mauro Maggioni
% mauro@math.duke.edu
% @ Copyright 2006-2010  Duke University
%

%
% Handle parameters
%
if isempty(cOpts), cOpts = []; end;

lEuclideanDistances = true;

if (~isfield(cOpts,'Dim')),
    if isempty(cX),
        cOpts.Dim = Inf;
    else
        cOpts.Dim = size(cX,1);
    end;
end;

if (~isfield(cOpts,'Sigma')) || (isempty(cOpts.Sigma)),                      cOpts.Sigma = cOpts.Delta; end;
if (~isfield(cOpts,'DownSample')) || (isempty(cOpts.DownSample)),            cOpts.DownSample = true; end;
if (~isfield(cOpts,'DoSVD'))  || (isempty(cOpts.DoSVD)),                     cOpts.DoSVD = false; end;
if (~isfield(cOpts,'DoD'))  || (isempty(cOpts.DoD)),                         cOpts.DoD = false; end;
if (~isfield(cOpts,'SVDDelta')) || (isempty(cOpts.SVDDelta)),                cOpts.SVDDelta = cOpts.Delta; end;
if (~isfield(cOpts,'NormalizedSVD')) || (isempty(cOpts.NormalizedSVD)),      cOpts.NormalizedSVD = 1; end;
if (~isfield(cOpts,'KeepV')) || (isempty(cOpts.KeepV)),                      cOpts.KeepV = 1; end;
if (~isfield(cOpts,'AvoidSmallClusters')) || (isempty(cOpts.AvoidSmallClusters)), cOpts.AvoidSmallClusters = 0; end;
if (~isfield(cOpts,'Randomized')) || (isempty(cOpts.Randomized)),            cOpts.Randomized = 1; end;
if (~isfield(cOpts,'DistInfo')) || (isempty(cOpts.DistInfo)),                cOpts.DistInfo = []; else lEuclideanDistances = false; end;
if (~isfield(cOpts,'SVDDim')) || (isempty(cOpts.SVDDim)),                    cOpts.SVDDim = cOpts.Dim;
elseif (cOpts.SVDDim > cOpts.Dim) && (lEuclideanDistances),              cOpts.SVDDim = cOpts.Dim; end;
if (~isfield(cOpts,'NNInfo')),                                               cOpts.NNInfo = []; end;
if (~isfield(cOpts,'RandomizedSVDThres')) || (isempty(cOpts.RandomizedSVDThres)), cOpts.RandomizedSVDThres = Inf; end;
if (~isfield(cOpts,'Verbose')), cOpts.Verbose = 1; end;
if (~isfield(cOpts,'FastNNSearcher')), cOpts.FastNNSearcher = []; end;
if (~isfield(cOpts,'KeepCenter')) || (isempty(cOpts.KeepCenter)),                      cOpts.KeepCenter = 0; end;
if (~isfield(cOpts,'RandomizedSVD')), cOpts.RandomizedSVD = false; end;
if (~isfield(cOpts,'kNN')), cOpts.kNN = []; end;    


% Get dimensionality and number of points
lDim = cOpts.Dim;
if isempty(cX),
    if isempty(cOpts.DistInfo),
        fprintf('\nFastRoughSetDownsample: Error: neither cX nor cOpts.DistInfo provided.\n');
        return;
    end;
    lNumberOfPoints = length(cOpts.DistInfo.idxs);
else
    lNumberOfPoints = size(cX,2);
end;

if (~isfield(cOpts,'RestrictSetIdxs')) || (isempty(cOpts.RestrictSetIdxs)),  cOpts.RestrictSetIdxs = uint32(1:lNumberOfPoints); end;

% Get rid of trivial case
if lNumberOfPoints==1,
    vNet.Idxs = uint32(1); vNet.InvIdxs = uint32(1); vNet.D = 0; vNet.S = zeros(1,cOpts.SVDDim);
    if cOpts.KeepV,
        if lEuclideanDistances,
            vNet.V = zeros([1,lDim,cOpts.SVDDim],'single');
        else
            vNet.V = zeros([1,cOpts.SVDDim,cOpts.SVDDim],'single');
        end;
    end;
    return;
end;

% Allocate memory
vNet.Idxs(lNumberOfPoints,1) = uint32(0);
vNet.InvIdxs(lNumberOfPoints,1) = uint32(0);

if (cOpts.DoSVD),
    % Allocate memory for SVD results
    vNet.S(lNumberOfPoints,cOpts.SVDDim) = double(0);
    if cOpts.KeepV
        if lEuclideanDistances,
            %vNet.V = sparse(lNumberOfPoints,lDim,cOpts.SVDDim);
            vNetVSize = [lDim,cOpts.SVDDim];
        else
            %vNet.V = sparse(lNumberOfPoints,cOpts.SVDDim,cOpts.SVDDim);
            vNetVSize = [cOpts.SVDDim,cOpts.SVDDim];
        end;
        vNet.V = {};
    end;
    if cOpts.KeepCenter && lEuclideanDistances, % added by YM
        vNet.Center(lNumberOfPoints,lDim) = double(0);
    end;
else
    % No svd results if not asked for them
    vNet.S = [];
    if cOpts.KeepV,
        vNet.V = [];
    end;
end;
if (cOpts.DoD), vNet.D(lNumberOfPoints,1) = double(0); else vNet.D = []; end;

lNotChosenIdxs = uint32(sort(cOpts.RestrictSetIdxs));

% Index of the first point to be chosen
lChosenIdx = uint32(1);
% Set up options for nearest neighbor search
lOptsForNrsearch.DistInfo       = cOpts.DistInfo;
lOptsForNrsearch.NNInfo         = cOpts.NNInfo;
lOptsForNrsearch.FastNNSearcher = cOpts.FastNNSearcher;
lOptsForNrsearch.SaveOpts       = true;
lOptsForNrsearch.XIsTransposed  = false;

lS = [];
lV = [];

lStatusString = [];

lRestrictDistInfoToSubSetOpts = struct('Idxs_sorted',true);

if ~lEuclideanDistances,
    global DistInfo;
    DistInfo = lOptsForNrsearch.DistInfo;
    if ~isfield(lOptsForNrsearch.DistInfo,'isinidxs'),
        lOptsForNrsearch.DistInfo = AddSortedIdxsToDistInfo( lOptsForNrsearch.DistInfo );
    end;
end;

if strcmpi(lOptsForNrsearch.FastNNSearcher,'nn_search'),
    cXt = cX';
    lAtria = nn_prepare( cXt );
end;

while ~isempty(lNotChosenIdxs),
    %% Select the next point
    if ~cOpts.Randomized,
        vNet.Idxs(lChosenIdx) = lNotChosenIdxs(1);
    else
        lRandIdx = floor(rand(1)*length(lNotChosenIdxs))+1;%   randi(length(lNotChosenIdxs),1);
        vNet.Idxs(lChosenIdx) = lNotChosenIdxs(lRandIdx);
    end;
    %% Find its cOpts.Delta-neighborhood or cOpts.kNN neighbors
    if ~strcmpi(lOptsForNrsearch.FastNNSearcher,'nn_search')
        if isempty(cOpts.kNN),
            [lNNCount, lNN, lNNDists, lNNInfo] = nrsearch( cX,vNet.Idxs(lChosenIdx),0,max([cOpts.Delta,cOpts.SVDDelta]),[],lOptsForNrsearch );
        else
            [lNNCount, lNN, lNNDists, lNNInfo] = nrsearch( cX,vNet.Idxs(lChosenIdx),cOpts.kNN,[],[],lOptsForNrsearch );
            cOpts.Delta = max(lNNDists{1});
            cOpts.SVDDelta = cOpts.Delta;
        end;
    else
        if isempty(cOpts.kNN),
            [lNNCount, lNeighbors] = range_search( cXt,lAtria,double(vNet.Idxs(lChosenIdx)),max([cOpts.Delta,cOpts.SVDDelta]), -1 );
            lNNInfo = []; clear lNNDists;
            for p = size(lNeighbors,1):-1:1,
                [lNNDists{p},lSortedIdxs] = sort(lNeighbors{p,2});
                lNN{p} = lNeighbors{p,1}(lSortedIdxs);
            end;
        else
            [lNeighbors, lDists] = nn_search( cXt,lAtria,double(vNet.Idxs(lChosenIdx)),cOpts.kNN, -1 );
            lNNInfo = []; clear lNNDists; 
            lNNDists{1} = lDists; 
            lNN{1} = lNeighbors;
            cOpts.Delta = max(lNNDists{1});
            cOpts.SVDDelta = cOpts.Delta;
        end;
    end;
    if lOptsForNrsearch.SaveOpts,
        lOptsForNrsearch.ReuseOpts = true;
        lOptsForNrsearch.SaveOpts = false;
    end;
    % Indices and distances of the points in the neighborhood
    lCloseIdxs = [];
    if iscell(lNNDists),
        if length(lNNDists)>0,
            lNNDists=lNNDists{1};
            lCloseIdxs = uint32(lNN{1}); %union(lNN{1},lNotChosenIdxs(1));
        end;
    else
        lCloseIdxs = uint32(lNN);
    end;
    % If SVD is going to be computed on larger ball, correct nearby points accordingly
    % SVD neighbors
    lCloseSVDIdxs = lCloseIdxs;
    if cOpts.SVDDelta>cOpts.Delta,
        % Net neighbors
        lIdxs = find(lNNDists<=cOpts.Delta);
        lCloseIdxs = lCloseIdxs(lIdxs);
        lNNDists = lNNDists(lIdxs);
    end;
    if isempty(lOptsForNrsearch.NNInfo),
        lOptsForNrsearch.NNInfo = lNNInfo;                    % May contain useful info at future iterations, such as Atria structure if TSTool is used.
    end;
    % Save local volume
    vNet.count(lChosenIdx) = length(lCloseIdxs);
    
    %% Compute the ``local degree'', if requested
    if (cOpts.DoD),
        % Compute local weights
        lCloseW = exp(-(lNNDists/cOpts.Sigma).^2);
        % Assign the degree to the vertex.
        vNet.D(lChosenIdx) = sum(lCloseW);
    end;
    
    %% Compute the local SVD if requested
    if cOpts.DoSVD,
        if cOpts.SVDDelta==cOpts.Delta,
            lCloseSVDIdxs = lCloseIdxs;
        end;
        vNet.SVDcount(lChosenIdx) = length(lCloseSVDIdxs);
        
        % Check if need to perform randomized SVD
        if length(lCloseSVDIdxs)>cOpts.RandomizedSVDThres,
            lTempIdxs = randperm(length(lCloseSVDIdxs));
            lCloseSVDIdxs = lCloseSVDIdxs(lTempIdxs(1:cOpts.RandomizedSVDThres)); clear lTempIdxs;
        end;
        
        if length(lCloseSVDIdxs)>1,
            if lEuclideanDistances,
                % Center the points
                lXlocal = cX(:,lCloseSVDIdxs);
                if ~cOpts.KeepCenter                    
                    lXlocal = bsxfun(@minus,lXlocal,mean(lXlocal,2));
                else
                    vNet.Center(lChosenIdx,:) = mean(lXlocal,2);                    
                    lXlocal = bsxfun(@minus,lXlocal,vNet.Center(lChosenIdx,:));
                end;                
                if cOpts.KeepV,
                    if isempty(lS)
                        [lV lS] = RandPCA(lXlocal, min([min(size(lXlocal)),cOpts.SVDDim]));
                    end;
                    lActualDim = min(cOpts.SVDDim,size(lV,2));
                    vNet.V{lChosenIdx} = zeros(vNetVSize,'single');
                    vNet.V{lChosenIdx}(1:size(lV,1),1:lActualDim) = single(lV(:,1:lActualDim));                         %vNet.V(lChosenIdx,1:size(lV,1),1:lActualDim) = sparse(lV(:,1:lActualDim));
                    if size(lS,2)>1, lS=diag(lS);  else  lS=lS(1,1);    end;
                else
                    if isempty(lS),
                        if lDim > cOpts.SVDDim
                            lS = svds(lXlocal, cOpts.SVDDim);
                        else
                            if ~cOpts.RandomizedSVD
                                lS = svd(lXlocal,0);
                            else
                                [lV,lS] = RandPCA(lXlocal,min([size(lXlocal),cOpts.SVDDim]));
                                lS=diag(lS);
                            end;
                        end
                    end;
                end;
            else
                % First of all construct the local dissimilarity matrix
                lDistLocal = FastDistRestrict( lOptsForNrsearch.DistInfo, lCloseSVDIdxs );
                % Do Classical MultiDimensionalScaling
                [lTmp,lEigInfo] = ClassicalMDS( lDistLocal.^2, struct('kEigenVecs',cOpts.SVDDim) );
                lS = lEigInfo.EigenVals;
                if cOpts.KeepV,
                    lV = lEigInfo.EigenVecs;
                end;
            end;
            if cOpts.NormalizedSVD,
                lS = lS./((sqrt(length(lCloseSVDIdxs))));
            end;
            if length(lS)>cOpts.SVDDim,
                lS = lS(1:cOpts.SVDDim);
            end;
            
            % Save the singular values
            vNet.S(lChosenIdx,1:length(lS)) = lS.^2;
        else
            % Save the singular values
            vNet.S(lChosenIdx,1:cOpts.SVDDim) = 0;
            if cOpts.KeepV,
                vNet.V{lChosenIdx} = zeros(vNetVSize,'single');
            end;
            % vNet.V is already taken care of
            if lEuclideanDistances && cOpts.KeepCenter
                % Center the points
                vNet.Center(lChosenIdx,:) = cX(:,lNotChosenIdxs(1));
            end
        end;
        
        lS = [];
        lV = [];
    end;
    
    %% Downsample
    if cOpts.DownSample,
        % The point selected and its cOpts.Delta-neighbors shouldn't be selected again: take them off the list of available good points
        lCloseIdxs = sort([lCloseIdxs,vNet.Idxs(lChosenIdx)]);
         
        lNotChosenIdxs = setdiff(lNotChosenIdxs,lCloseIdxs,'sorted');
         
        %BBCREVISIST Changed the line - setdiff_sorted was not in the code
         %from Mauro
        % lNotChosenIdxs = setdiff_sorted(lNotChosenIdxs,lCloseIdxs);
      
        %lNotChosenIdxs = uint32(lNotChosenIdxs(find(~(ismember_fast_sorted(lNotChosenIdxs,union(lCloseIdxs,lNotChosenIdxs(1)))))));
    else
        lNotChosenIdxs(lRandIdx)=[];
    end;
    
    %% Verbose
    if cOpts.Verbose,
        if rem(lChosenIdx,50)==0,
            fprintf(repmat(['\b'],[1,length(lStatusString)]));
            lStatusString = sprintf('%d net pts., %d total remaining (out of %d)',lChosenIdx, length(lNotChosenIdxs), lNumberOfPoints);
            fprintf(lStatusString);
        end;
    end;
    
    % Save the delta at this scale
    if ~isempty(cOpts.kNN),
        vNet.Delta(lChosenIdx) = cOpts.Delta;
    end;
    
    % Go to the next point
    lChosenIdx = lChosenIdx+1;                                              % For debug purpose only: figure(1);plot(length(vNet.Idxs),length(lNotChosenIdxs),'o');hold on;drawnow;pause
end;

if isempty(cOpts.kNN),
    vNets.Delta = cOpts.Delta;
end;

if cOpts.Verbose,
    fprintf(repmat(['\b'],[1,length(lStatusString)]));
    lStatusString = sprintf('%d net pts., %d total remaining (out of %d)',lChosenIdx-1, length(lNotChosenIdxs), lNumberOfPoints);
    fprintf(lStatusString);
end;

vNet.NNInfo = lNNInfo;

%% Some post-processing steps
% Get rid of part of pre-allocated memory never used.
vNet.Idxs(lChosenIdx:length(vNet.Idxs))=[];


%% Find the "inverse" map
if isempty(lOptsForNrsearch.DistInfo),
    [lCount,vNet.InvIdxs] = nrsearch( cX(:,vNet.Idxs), cX, 1, [], [], struct('ReturnAsArrays',1,'XIsTransposed',lOptsForNrsearch.XIsTransposed,'ReuseOpts',false,'SaveOpts',false,'FastNNSearcher',lOptsForNrsearch.FastNNSearcher) );
else
    lOptsForNrsearch.ReuseOpts = false;
    lOptsForNrsearch.SaveOpts = false;
    vNet.InvIdxs = FastAssignToCenters( vNet.Idxs, uint32(1:lNumberOfPoints), lOptsForNrsearch.DistInfo );
end;

%% Clean up
if (cOpts.DoD),
    vNet.D(lChosenIdx:length(vNet.InvIdxs))=[];
end;
if cOpts.DoSVD,
    vNet.S(lChosenIdx:length(vNet.InvIdxs),:) = [];
    if cOpts.KeepV,
        lV = vNet.V;
        vNet.V = zeros([length(lV),vNetVSize],'single');
        for k = 1:length(lV),
            if ~isempty(lV{k}),
                vNet.V(k,:,:) = single(lV{k});
            end;
        end;
        clear lV;
        %vNet.V(lChosenIdx:length(vNet.InvIdxs),:,:) = [];
    end;
end;
if cOpts.KeepCenter
    vNet.Center(lChosenIdx:end,:) = [];
end

%% More post-processing: removes small clusters if desired
if cOpts.AvoidSmallClusters>0,
    lMedianCount = median(vNet.count);
    % Get the distribution of the sizes of the clusters
    [lSortedCounts,lSortedIdxs] = sort( vNet.count );
    % Find the small clusters
    lSmallClusterIdxs = max(find(lSortedCounts<cOpts.AvoidSmallClusters * lMedianCount));
    if length(lSmallClusterIdxs)>0,
        % Discard the small clusters
        idxToBeRemoved = lSortedIdxs(1:lSmallClusterIdxs);
        vNet.Idxs(idxToBeRemoved) = [];
        if cOpts.DoSVD,
            vNet.S(idxToBeRemoved,:) = [];
            if cOpts.KeepV,
                vNet.V(idxToBeRemoved,:,:) = [];
            end;
            vNet.SVDcount(idxToBeRemoved)=[];
        end;
        if (cOpts.DoD),
            vNet.D(idxToBeRemoved) = [];
        end;
        vNet.count(idxToBeRemoved) = [];
        
        % Update the "inverse" map
        [lTmpCount, vNet.InvIdxs] = nrsearch( cX(:,vNet.Idxs), cX,1, [], [], struct('ReturnAsArrays',1) );
        
        % Update the count of points associated with each center
        vNet.count = hist(vNet.InvIdxs,1:length(vNet.Idxs));
    end;
end;

if cOpts.Verbose,
    fprintf('\n');
end;

return;
