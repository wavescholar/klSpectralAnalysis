function vNets = FastRoughSetMultiscaleSVD( cX, cOpts )

%
% function vNets = FastRoughSetMultiscaleSVD( cX, cOpts )
%
% IN:
%   cX      : D by N matrix of N points in D dimensions?
%   cOpts   : structure of options, with the following fields:
%               Epsilon        : radius for the finest scale.
%               Delta          : if scalar, dilation factor to go from scale to scale, otherwise vector of scales.
%               SVDDelta       : if scalar, multiplies Delta, otherwise vector of scales. Default: 1.
%               kNN            : if scalar, number of nn's added from scale to scale, otherwise vector of nn's to be used at all scales. Overwrites Delta and SVDDelta. Default: [].
%               NumberOfScales : number of scales
%               ApproxSVD      : 0,1
%               NormalizedSVD  : 0,1
%               [DecreasingNets] : 0,1: whether the net at scale j is contained in the one at scale j-1 or not. Default: 1.
%               [DistInfo]     : passed to FastRoughSetDownsample
%               [MinNetPts]    : stop if net has less than this number of points. Default: 0.
%               [MinNetDeltaPts]: stop if net at current scale j has less than this number of points less than the net
%                                at scale j-5. Default: 0.
%               [MinNetDeltaLogRadius] : stop if net at current scale j has a radius r_j s.t. r_{j-1}/r_j>MinNetDeltaLogRadius
%             
%             NOTE: Any other members of this structure are passed to FastRoughSetDownsample
%
% OUT:
%   vNets   : structure array indexed by scale (so it's cOpts.NumberOfScales long), with the following fields:
%               AbsIdxs        : indices into cX of the points of the net.
%               AbsInvIdxs     : M vector s.t. cX(InvIdxs,:) is set-wise equal to the net. In other words this maps cX onto the net.
%               Idxs           : indices into net at the previous scale (j-1), if cOpts.DecreasingNets==1, otherwise equal to AbsIdxs.
%               InvIdxs        : same as InvIdxs, but relative to the previous scale, if cOpts.DecreasingNets==1, otherwise equal to AbsInvIdxs.
%               S              : |Idxs| by N matrix of singular values at each point in the net
%               V              : |Idxs| by N by N tensor of singular vectors at each point in the net
%               D              : |Idxs| vector of weights of each point of the net
%               count          : number of neighbors used for local svd computation
%               Delta          : vector of scales
%
% EXAMPLE:
%  lN = 1000; lParam=linspace(0,2*pi,lN);cX=[sin(lParam)',cos(lParam)']+0.01*randn(lN,2);cX(1:1000,:)=cX(1:1000,:)+0.1*rand(1000,2);
%  lOpts = struct('Epsilon',0.05,'Delta',2,'SVDDelta',1.5,'NumberOfScales',5,'ApproxSVD',0);
%  [DistInfo.count,DistInfo.idxs, DistInfo.dists] = nrsearch(cX', cX', 100);
%  lOpts = struct('Epsilon',0.05,'Delta',2,'SVDDelta',1.5,'NumberOfScales',5,'ApproxSVD',0,'DistInfo',DistInfo);
%  tic;vNets = FastRoughSetMultiscaleSVD( cX, lOpts );toc;
%  for lj = 1:lOpts.NumberOfScales; DisplayNetWithSVD( cX, vNets(lj).AbsIdxs, vNets(lj).S, vNets(lj).V );end;
%
%
% USES:
%   FastRoughSetDownsample
%
% Mauro Maggioni
% mauro@math.duke.edu
% @ Copyright 2006  Duke University
%


%% Parameter check & initialization
if ~isfield(cOpts,'DecreasingNets'),            cOpts.DecreasingNets = 1;       end;
if ~isfield(cOpts,'SVDDelta'),                  cOpts.SVDDelta = 1;             end;
if ~isfield(cOpts,'kNN'),                       cOpts.kNN = [];                 end;
if ~isfield(cOpts,'XIsTransposed'),             cOpts.XIsTransposed = false;    end;
if ~isfield(cOpts,'MinNetPts'),                 cOpts.MinNetPts = 0;            end;
if ~isfield(cOpts,'MinNetDeltaPts'),            cOpts.MinNetDeltaPts = 0;       end;
if ~isfield(cOpts,'MinNetDeltaLogRadius'),      cOpts.MinNetDeltaLogRadius = 1; end;
if ~isfield(cOpts,'Verbose'),                   cOpts.Verbose = 1;              end;

[D,N] = size(cX);

lkNN = [];

if isempty(cOpts.kNN),
    if (length(cOpts.Delta)==1),
        lDeltas = zeros(cOpts.NumberOfScales,1);
        for lj = 1:cOpts.NumberOfScales,
            lDeltas(lj) = (cOpts.Delta.^(lj-1))*(cOpts.Epsilon);
        end;
    else
        lDeltas = cOpts.Delta*cOpts.Epsilon;
    end;
    
    if (length(cOpts.SVDDelta)==1) || (length(cOpts.SVDDelta)~=cOpts.NumberOfScales),
        lSVDDeltas = cOpts.SVDDelta * lDeltas;
    else
        lSVDDeltas = cOpts.SVDDelta;
    end;
else
    if (length(cOpts.kNN)==1),
        lkNN = zeros(cOpts.NumberOfScales,1);
        for lj = 1:cOpts.NumberOfScales,
            lkNN(lj) = min([lj*cOpts.kNN,N]);            
        end;
    else
        lkNN = cOpts.kNN;
    end;
end;


%% Build the first scale (finest)
if cOpts.Verbose,
    if isempty(cOpts.kNN),
        fprintf('Scale 1/%d, radius %.4f...',cOpts.NumberOfScales,lDeltas(1));
    else
        fprintf('Scale 1/%d, kNN %d...',cOpts.NumberOfScales,lkNN(1));
    end;
end;

lOpts = cOpts;
lOpts.DoSVD = true;
if isempty(lkNN),
    lOpts.Delta = lDeltas(1);
    lOpts.SVDDelta = lSVDDeltas(1);
else
    lOpts.kNN = lkNN(1);
end;
vNets{cOpts.NumberOfScales} = [];
vNets{1} = FastRoughSetDownsample( cX, lOpts );
vNets{1}.AbsIdxs = vNets{1}.Idxs;
vNets{1}.AbsInvIdxs = vNets{1}.InvIdxs;
vNets{1}.Delta = lOpts.Delta;
lOpts.NNInfo = vNets{1}.NNInfo;                         % Re-use fast nearest neighbor structure if possible

%% Check if approximate or exact svd is requested
if cOpts.ApproxSVD,
    %% Approximate svd, computed in multiscale fashion
    for lj = 2:cOpts.NumberOfScales, 
        if isempty(lkNN),
            lOpts.Delta = lDeltas(lj);
            lOpts.SVDDelta = lSVDDeltas(lj);
        else
            lOpts.kNN = lkNN(lj);
        end;
        if cOpts.Verbose,
            if isempty(cOpts.kNN),
                fprintf('Scale %d/%d, radius=%.4f...',lj,cOpts.NumberOfScales, lOpts.Delta);
            else
                fprintf('Scale %d/%d, kNN=%d...',lj,cOpts.NumberOfScales, lOpts.kNN);
            end;
        end;
        if cOpts.XIsTransposed,
            vNets{lj} = FastRoughSetDownsample( cX(vNets{lj-1}.AbsIdxs,:), lOpts ); 
        else
            vNets{lj} = FastRoughSetDownsample(cX(:, vNets{lj-1}.AbsIdxs), lOpts ); 
        end;      
        vNets{lj}.AbsIdxs = vNets{lj-1}.AbsIdxs(vNets{lj}.Idxs);
        vNets{lj}.AbsInvIdxs = vNets{lj}.InvIdxs(vNets{lj-1}.AbsInvIdxs);
        vNets{lj}.Delta = lOpts.Delta;
        % If the number of points in the net falls below threshold, or the running average size of the nets is not decrease
        if (length(vNets{lj}.count) < cOpts.MinNetPts) || ((lj>5) && ((((length(vNets{lj-5}.count)-length(vNets{lj}.count))<cOpts.MinNetDeltaPts) && (cOpts.MinNetDeltaPts>0)) || (lDeltas(lj-1)/lDeltas(lj)>cOpts.MinNetDeltaLogRadius))),
            break;
        end; 
    end;
else
    %% True svd
    for lj = 2:cOpts.NumberOfScales,
        if isempty(lkNN),
            lOpts.Delta = lDeltas(lj);
            lOpts.SVDDelta = lSVDDeltas(lj);
        else
            lOpts.kNN = lkNN(lj);
        end;
        if cOpts.Verbose,
            if isempty(cOpts.kNN),
                fprintf('Scale %d/%d, radius=%.4f...',lj,cOpts.NumberOfScales, lOpts.Delta);
            else
                fprintf('Scale %d/%d, radius=%.4f...',lj,cOpts.NumberOfScales, lOpts.kNN);
            end;
        end;
        if cOpts.DecreasingNets,
            lOpts.RestrictSetIdxs = vNets{lj-1}.AbsIdxs;
        end;
        vNets{lj} = FastRoughSetDownsample( cX, lOpts );
        vNets{lj}.AbsIdxs = vNets{lj}.Idxs;       
        vNets{lj}.AbsInvIdxs = vNets{lj}.InvIdxs;
        vNets{lj}.Delta = lOpts.Delta;
        % If the number of points in the net falls below threshold, or the running average size of the nets is not decrease
        if (length(vNets{lj}.count) < cOpts.MinNetPts) || ((lj>5) && ((((length(vNets{lj-5}.count)-length(vNets{lj}.count))<cOpts.MinNetDeltaPts) && (cOpts.MinNetDeltaPts>0)) || ((isempty(lkNN)) && (lDeltas(lj-1)/lDeltas(lj)>cOpts.MinNetDeltaLogRadius)))),
            break;
        end;        
    end;
end;

vNets(lj:cOpts.NumberOfScales)=[];

vNets = cell2mat(vNets);

return;