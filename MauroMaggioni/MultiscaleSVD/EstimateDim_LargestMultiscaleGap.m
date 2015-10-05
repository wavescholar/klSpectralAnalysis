function vStats = EstimateDim_LargestMultiscaleGap( cX, cData, cOpts )

%
% function vStats = EstimateDim_LargestMultiscaleGap( cX, cData, cOpts )
%
%
% IN:
%   cX          : D by N data matrix of N points in R^D
%   cData       : structure containing data information:
%                   Stats   :   structure containing the following fields:
%                       S       :   (#scales)*(#points)*(#dimensions) tensor whose (i,j,k) entry is the k-th singular value of PCA at scale i for point j
%                   Nets    :   (#scales) array of structures with the following fields:
%                       count   :   (#net points) vector with the number of points in each cell at scale j
%                   NetsOpts:   structure containing the following fields:
%                       Delta   : (#scales) vector containing the set of scales
%   [cOpts]     : structure containing the following fields:
%                   [Verbose]       : Default: false.
%                   [Ptwise]        : Pointwise vs. global dimensionality. Default: false.
%                   [UseSmoothedS]  : Uses smoothed multiscale singular values. Default: false.
%
% OUT:
%   vStats      : structure with the following fields:
%                   GoodScales  : if cOpts.Ptwise==false, a vector of indices (into NetsOpts.Delta) of good scales, 
%                                 if cOpts.Ptwse==true,   a N cell array where the i-th cell contains a vector of indices (into NetsOpts.Delta) of good scales
%                   DimEst      : estimated dimensionality: a real number (if ~cOpts.Ptwise) or a N vector (if cOpts.Ptwise)
%                   
%
% EXAMPLE
%   vStats = EstimateDim_LargestMultiscaleGap( X, Data );
%   j=7;figure;PlotWithStd(1:size(Data.Nets(1).S,2)-1,squeeze(vStats.Gaps(j,Data.Nets(j).AbsInvIdxs,:))');
%   j=5;figure;PlotWithStd(1:size(Data.Nets(1).S,2),squeeze(Data.Stats.Tychonov.S(j,Data.Nets(j).AbsInvIdxs,:))');hold on;PlotWithStd(1:size(Data.Nets(1).S,2),squeeze(Data.Nets(j).S(Data.Nets(j).AbsInvIdxs,:))','r-');axis tight;
%
%
%

% (c) Copyright 
% Mauro Maggioni
% Duke University, 2009
%

if nargin<3,                        cOpts = [];                 end;
if ~isfield(cOpts,'Verbose'),       cOpts.Verbose = false;      end;
if ~isfield(cOpts,'Ptwise'),        cOpts.Ptwise  = false;      end;
if ~isfield(cOpts,'UseSmoothedS'),  cOpts.UseSmoothedS = false; end;

vStats.DimEst = 0;

% Number of points
lN = size(cX,2);

for j = length(cData.Nets):-1:1,
    lVolumes(j) = mean(cData.Nets(j).count);
end;

if ~cOpts.Ptwise
    % Compute the mean gaps between the singular values as a function of scale
    if cOpts.UseSmoothedS,
        if isfield(cData.Stats,'Tychonov'),
            lMeanS = squeeze(mean(cData.Stats.Tychonov.S,2));
        else
            warning('EstimateDim_LargeMultiscaleGap: UseSmoothedS option selected, but no smoothed singular values computed.');
            lMeanS = squeeze(mean(cData.Stats.S,2));
        end;
    else
        lMeanS = squeeze(mean(cData.Stats.S,2));
    end;
    % For debug purposes: figure;for k = 1:size(cData.Stats.S,3);PlotWithStd(lDeltas,squeeze(cData.Stats.S(:,:,k)));hold on;end;axis tight;figure;for k = 1:size(cData.Stats.S,3);PlotWithStd(lMeanS(:,1),squeeze(cData.Stats.S(:,:,k)));hold on;end;axis tight;
    vStats      = EstimateDimFromSpectra( cData.NetsOpts.Delta', lVolumes, cOpts, lMeanS );
else
    GoodScales{lN} = [];
    DimEst(lN)     = 0;
    for k = 1:lN,
        if cOpts.UseSmoothedS,
            if isfield(cData.Stats,'Tychonov'),
                lS = squeeze(cData.Stats.Tychonov.S(:,k,:));
            else
                warning('EstimateDim_LargeMultiscaleGap: UseSmoothedS option selected, but no smoothed singular values computed.');
                lS = squeeze(cData.Stats.S(:,k,:));
            end;
        else
            lS = squeeze(cData.Stats.S(:,k,:));
        end;
        lStats        = EstimateDimFromSpectra( cData.NetsOpts.Delta', lVolumes, cOpts, lS );
        GoodScales{k} = lStats.GoodScales;
        DimEst(k)     = lStats.DimEst;        
    end;
    vStats.GoodScales = GoodScales;
    vStats.DimEst     = DimEst;
end;

return;