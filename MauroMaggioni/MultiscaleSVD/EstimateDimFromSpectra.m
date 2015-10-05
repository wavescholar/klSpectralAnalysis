function vStats = EstimateDimFromSpectra( cDeltas, cVolumes, cOpts, S_MSVD )

%%
%
% Estimates intrinsic dimensionality and good scales given the multiscale singular values
%
% IN:
%   cDeltas         : (#scales) vector of scales
%   cVolumes        : (#scales) vector of average volume of cells at each scale
%   [cOpts]         : structure containing the following fields:
%                       [Verbose]       : Default: false.
%   S_MSVD          : (#scales)*(#dimensions) matrix of singular values: the (i,j) entry is the j-th singular value of a cell at scale cDeltas(i) around a point
%
% OUT:
%   vStats          : structure containing the following fields:
%                       DimEst      : estimate of intrinsic dimensionality
%                       GoodScales  : vector of indices (into the vector cDeltas) of good scales used for dimension estimate
%
%
% (c) Copyright
% Mauro Maggioni, Duke University, 2008
%
%%

warning off;

if ~isfield(cOpts,'Verbose'),   cOpts.Verbose = false; end;

% Parameters
epsnoise = 1e-3;            % Parameter for determining possible noise
pBacktrackDim = 2;          % How many dimensions to backtrack at once when estimating the dimension
pMinNeighbors = 10;         % Don't look at scales that contain less than this number of neighbors

lMinScaleByNeighbors = find(cVolumes>=pMinNeighbors);

% Deltas
lDeltas     = cDeltas(lMinScaleByNeighbors);
lVolumes    = cVolumes(lMinScaleByNeighbors);
lMeanS      = S_MSVD(lMinScaleByNeighbors,:);
lMeanGapS   = [abs(diff(lMeanS,1,2)),lMeanS(:,max([1,size(lMeanS,2)-1]))];

% Number of scales
lJ = size(lMeanS,1);
% Number of singular values
lS = size(lMeanS,2);

%% Start by estimating the dimension from large scales (this should be an upper bound of the actual dimension)
lLargeScales = find(lDeltas>0.85*max(lDeltas));%find(lDeltas>0.75*max(lDeltas));
lMaxDim = EstimateMaxDim( lDeltas, lLargeScales,lMeanS,lMeanGapS );
if cOpts.Verbose, fprintf('\n\t\t Max dim: %d...',lMaxDim); end;

%% Compute the mean of the noise multiscale spectra
lMeanNoiseSpectra   = mean(lMeanS(:,min([lMaxDim+1,lS]):lS),2);
lMaxNoiseSpectra    = lMeanS(:,min([lMaxDim+1,lS]));
if lMaxDim>=lS,
    lMeanNoiseSpectra = 0*lMeanNoiseSpectra;
    lMaxNoiseSpectra  = 0*lMaxNoiseSpectra;
end;
lNoiseSigma         = mean(lMeanNoiseSpectra(max(lJ-5,1):lJ));
if lNoiseSigma>epsnoise*max(lMeanS(:,lMaxDim)),
    lIdxNoiseDecreasing = find(abs(lMeanNoiseSpectra-lNoiseSigma)/lNoiseSigma<0.2);
    if lMaxDim<size(lMeanS,2),
        lIdxNoiseIsSmall    = find(lMeanS(:,lMaxDim)-lMeanS(:,lMaxDim+1)>0.2*lMeanS(:,lMaxDim));
    else
        lIdxNoiseIsSmall = 1:lJ;
    end;
else
    lIdxNoiseDecreasing = 1:lJ;
    lIdxNoiseIsSmall = 1:lJ;
end;
%% Good scales for noise, both small and decreasing
lIdxGoodNoiseScales = lIdxNoiseIsSmall;
if (length(lIdxGoodNoiseScales)<5) && (~isempty(lIdxNoiseDecreasing)),
    lIdxGoodNoiseScales = lIdxNoiseDecreasing;
end;

%% Scales at which the manifold singular values are roughly constant
lMaxDimAsympt = mean(lMeanS(max(lJ-5,1):lJ,lMaxDim));
if lMaxDimAsympt>epsnoise*max(lMeanS(:,lMaxDim)),
    lIdxMaxDimNonAsympt = 1:max(find(abs(lMeanS(:,lMaxDim)-lMaxDimAsympt)/lMaxDimAsympt>epsnoise));
else
    lIdxMaxDimNonAsympt = 1:lJ;
end;

%% Scales at which the first singular value is growing
try
    tmp = lMeanS(lIdxGoodNoiseScales,1);
    lIdxGoodNoiseScales = lIdxGoodNoiseScales(find(abs(tmp-mean(tmp))<4*std(tmp)));
    [lFit,lFitStats]= polyfit(lDeltas(lIdxGoodNoiseScales),lMeanS(lIdxGoodNoiseScales,1),3);
    lFitVar = norm(lFitStats.R*pinv(lFitStats.R)')*lFitStats.normr/lFitStats.df;
    ldFit = diff(polyval(lFit,lDeltas(lIdxGoodNoiseScales)));
    lIncreasingTopSingVal = FindLargestContiguousConstantPiece( ldFit > -median(abs(ldFit)) ); % min([0,median(ldFit)]) );
catch
    lIncreasingTopSingVal = 1:length(lIdxGoodNoiseScales);
end;

%% Scales before jumps in the spectrum: look for jumps in the lMaxDim singular value
lDiffMeanS = diff(lMeanS(lIdxGoodNoiseScales,1:lMaxDim),1);
if size(lDiffMeanS,1)==1,
    lVarS = lDiffMeanS(:,1:lMaxDim);
else
    lVarS = var(lDiffMeanS(:,1:lMaxDim));
end;
lMinLargeJump = [];%size(lDiffMeanS,1);
for k = lMaxDim:lMaxDim,
    lMinLargeJump = lIdxGoodNoiseScales(max([lMinLargeJump,max(find(abs(lDiffMeanS(:,k)-mean(lDiffMeanS(:,k)))>100*lVarS(k)))]));
end;

vStats.GoodScales = intersect(intersect(lIdxGoodNoiseScales,lIdxMaxDimNonAsympt),lIdxGoodNoiseScales(lIncreasingTopSingVal));
if isempty(vStats.GoodScales),
    vStats.GoodScales = lIdxGoodNoiseScales;
end;
if ~isempty(lMinLargeJump),
    vStats.GoodScales = intersect(vStats.GoodScales,1:lMinLargeJump);
end;

if cOpts.Verbose, fprintf('\n\t\t Good scales: [%.3f,%.3f]...',min(lDeltas(vStats.GoodScales)),max(lDeltas(vStats.GoodScales)));end;

%% Finally, the smallest scale r_\min should be such that in average a ball of radius r_\min contains at least pi*lMaxDim points
lMinNumberOfPointsFound = false;
for j = 1:length(vStats.GoodScales),
    if mean(lVolumes(vStats.GoodScales(j)))>2*pi*lMaxDim,
        lMinNumberOfPointsFound = true;
        break;
    end;
end;

%% Find good scales and detect curvature
lGoodScalesFound = true;
vStats.GoodScales = vStats.GoodScales(j:length(vStats.GoodScales));
if cOpts.Verbose, fprintf('\n\t\t Good scales by considering sampling: [%.3f,%.3f]...',min(lDeltas(vStats.GoodScales)),max(lDeltas(vStats.GoodScales)));end;

vStats.DimEst = EstimateMaxDim( lDeltas, vStats.GoodScales,lMeanS(:,1:lMaxDim),lMeanGapS(:,1:lMaxDim) );

% Check linear vs. quadratic growth
if vStats.DimEst<size(lMeanS,2),
    lShiftedMeanS = lMeanS-repmat(lMeanS(:,vStats.DimEst+1),[1,size(lMeanS,2)]);
else
    lShiftedMeanS = lMeanS;
end;
lShiftedMeanGapS = [abs(diff(lShiftedMeanS,1,2)),lShiftedMeanS(:,max([1,size(lShiftedMeanS,2)]))]; %lMeanGapS;           %lShiftedMeanGapS(:,vStats.DimEst) = lMeanGapS(:,vStats.DimEst) - min(lMeanGapS(vStats.GoodScales,min([vStats.DimEst,size(lMeanGapS,2)])));

lEstDimChanged = true;
while (lEstDimChanged) && (vStats.DimEst >1),
    lEstDimChanged = false;
    if cOpts.Verbose, fprintf('\n\t\t Max dim.: %d',vStats.DimEst); end;
    % Check if (lMaxDim-1)-st gap is, at any scale, the largest gap
    lMaxGapsIdxs=almostmax(lShiftedMeanGapS(vStats.GoodScales,:),0.5);
    
    % Find the dominant gaps
    clear lSmallerDimIdxs;
    lSmallerDimIdxs{min([pBacktrackDim,vStats.DimEst-1])} = [];lCurDimIdxs = [];
    for i = 1:min([pBacktrackDim,vStats.DimEst-1]),
        for k = pBacktrackDim:-1:1,
            lSmallerDimIdxs{i}     = union(lSmallerDimIdxs{i},FindLargestContiguousConstantPiece(((lMaxGapsIdxs(:,k)==vStats.DimEst-i))));
            lCurDimIdxs            = union(lCurDimIdxs,FindLargestContiguousConstantPiece(((lMaxGapsIdxs(:,k)==vStats.DimEst))));
        end;
    end;
    
    if length(lCurDimIdxs)>3,
        for i = 1:min([pBacktrackDim,vStats.DimEst-1]),
            if (length(lSmallerDimIdxs{i})>3),
                lScales     = union(vStats.GoodScales(lSmallerDimIdxs{i}),vStats.GoodScales(lCurDimIdxs));
                lQuadFit    = polyfit(lDeltas(lScales),lShiftedMeanS(lScales,vStats.DimEst),2);
                if (lQuadFit(1)>0) && sum(lMaxGapsIdxs(1:floor(size(lMaxGapsIdxs,1)/2),1)==vStats.DimEst-i)>0.25*size(lMaxGapsIdxs,1),
                    vStats.DimEst  = vStats.DimEst-i;
                    lEstDimChanged = true;
                    break;
                end;
            end;
        end;
    else
        for i = 1:min([pBacktrackDim,vStats.DimEst-1]),
            if (length(lSmallerDimIdxs{i})>3)
                vStats.DimEst = vStats.DimEst-1;
                lEstDimChanged = true;
                break;
            end;
        end;
    end;
end;

% Update good range of scales
lMaxGapsIdxs      = almostmax(lShiftedMeanGapS,0.75);
lCorrectDimScales = find(sum(lMaxGapsIdxs(:,1:2)==vStats.DimEst,2)>0);
if vStats.DimEst<size(lShiftedMeanGapS,2),
    lDimGap           = (lShiftedMeanGapS(lCorrectDimScales,vStats.DimEst)-lShiftedMeanGapS(lCorrectDimScales,vStats.DimEst+1))./(lShiftedMeanGapS(lCorrectDimScales,vStats.DimEst));
    vStats.GoodScales = lCorrectDimScales(find(lDimGap>=median(lDimGap(find(lDimGap<=median(lDimGap))))));
end;
lContigIdxs       = FindLargestContiguousConstantPiece(diff(vStats.GoodScales));
vStats.GoodScales = vStats.GoodScales(min(lContigIdxs):max(lContigIdxs))+min(lMinScaleByNeighbors)-1;

if cOpts.Verbose,
    fprintf('\n\t\t Good scales: [%.3f,%.3f]...',min(lDeltas(vStats.GoodScales)),max(lDeltas(vStats.GoodScales)));
    fprintf('\n DimEst: %f',vStats.DimEst);
    fprintf('\n');
end;

if isempty(vStats.GoodScales),
    vStats.DimEst = vStats.DimEst+0.01;
    vStats.GoodScales = 1;
end;

warning on;

return;



%
% Find the largest interval where a function is nonzero
%
function vIdxs = FindLargestContiguousConstantPiece( cF )

vIdxs = [];
lLargestPieceLen = 0;

for k = 1:length(cF),
    lCurPieceLen = 0;
    if cF(k)~=0,
        for i = k+1:length(cF),
            if (cF(i)~=cF(k))
                if ~( (i<length(cF)-1) && (cF(i+1)==cF(k)) )
                    break;
                end;
            else
                lCurPieceLen = lCurPieceLen + 1;
            end;
        end;
        if lCurPieceLen > lLargestPieceLen,
            vIdxs = k:i-1;
            lLargestPieceLen = lCurPieceLen;
        end;
    end;
end;

return;




%
% Estimates upper bound on intrinsic dimension based on large scale singular values
%
function lMaxDim = EstimateMaxDim( lDeltas, lLargeScales,lMeanS,lMeanGapS )

epsnoise = 1e-3;

lLargeGapIdxs   = almostmax( lMeanGapS(lLargeScales,:), 0.5 );
lLargeGapIdxs   = unique( lLargeGapIdxs(:,1:min([size(lLargeGapIdxs,2),3])) );
lLargeGapSorted = sort(lLargeGapIdxs,'descend');
lRoughYScale    = max(max(lMeanS));
lRoughXScale    = max(lDeltas);
lMaxDim         = 0;

% Look for a linearly increasing singular value with a large gap, across large scales
for k = 1:length(lLargeGapSorted),
    lLinFit = polyfit(lDeltas(lLargeScales),lMeanS(lLargeScales,lLargeGapSorted(k)),1);
    if lLinFit(1)>epsnoise*lRoughYScale/lRoughXScale,
        lMaxDim = lLargeGapSorted(k);
        break;
    end;
end;

lMaxDim = max([max([max(max(lLargeGapSorted)),min([lMaxDim,size(lMeanS,2)-1])]),1]);

return;


function almostmaxidxs = almostmax( v, thres )

lSortedGaps = zeros(size(v));

% Sort the gaps at each scale
for k = 1:size(v,1),
    [lSortedValues,lSortedGaps(k,:)] = sort(v(k,:),'descend');
    lSmallIdxs = find(lSortedValues<=thres*lSortedValues(1));
    lSortedGaps(k,lSmallIdxs) = lSortedGaps(k,1);
end;

almostmaxidxs = lSortedGaps;

return;


% Find the largest gap as a function of scale
[lMax,lMaxIdxs] = max(v,[],2);
% Find the largest gap with the minimum index
lMinMaxIdx      = min(lMaxIdxs);

lvcomp = v(:,lMinMaxIdx);

for k = size(v,1):-1:1,
    lvrow = v(k,:);
    almostmaxidxs(k) = max(find(lvrow>=thres*lvcomp(k)));
end;

return;

function w = width( I )

w = max(I)-min(I);

return;