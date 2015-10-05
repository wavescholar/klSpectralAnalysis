function [EstDim,EstDimStats,Stats] = EstDim_MSVD( X, cOpts )

%
% function [EstDim,EstDimStats,Stats] = EstDim_MSVD( X, cOpts )
%
% IN:
%   X       : D by N matrix of N data points in D dimensions
%   [Opts]  : structure of options:
%               [MAXDIM]            : maximum intrinsic dimension of the data expected. Default: D.
%               [MAXAMBDIM]         : maximum ambient dimension of the data expected. Default: D.
%               [verbose]           : Default: false. Can be 1 or 2 for increasing levels of verbosity.
%               [Ptwise]            : estimates pointwise dimensionality. Default: false.
%               [PtIdxs]            : if Ptwise==true, this specifies at which points the intrinsic dimension should be estimated. Default: 0, meaning all points.
%               [NumberOfTrials]    : number of trials. These are run in parallel if possible. Default: 5.
%               [NetsOpts]          : to be used for constructing the nets. Default: computed automatically inside EstDim_MSVD.
%               [Deltas]            : will use this vector for the vector of scales. Default: computed adaptively inside EstDim_MSVD.
%               [kNNs]              : will use this vector for the vector of kNN's at different scales. Default: [].
%               [RandomizedSVDThres]: no more than this number of points will be used to estimate cov(X_{z,r}). Default: inf.
%               [UseSmoothedS]      : smooth or not the singular values. Default: false.
%               [KeepV]             : will keep the tangent vectors. Default: false
%               [AllPtsFast]        : will use external C code for computing the s.v.'s at all points [beta version]. Default: false.
%               [DownSample]        : will downsample points at each scale. Default: true.
%
%
% OUT:
%   EstDim      : estimated intrinsic dimension of the data: either a number of a vector of numbers if pointwise estimates are requested
%   EstDimStats : structure containing statistics computed during the estimation of intrinsic dimensionality
%                   DimEst      : vector of estimated intrinsic dimension during various randomized trials
%                   GoodScales  : an estimate of the range of good scales
%                   Timings: structure containing the following fields:
%                       EstDim_MSVD     : total running time
%                       MultiscaleNet   : time for computing the multiscale nets
%                       MultiscaleStats : time for computing the multiscale statistics
%                       EstimateDim     : time for estimating the intrinsic dimension
%   Stats        : structure containing the multiscale nets, SVDs, and multiscale statistics, as returned by FastRoughSetMultiscaleSVD
%                       
%
% USES:
%   Autoscales, FastRoughSetMultiscaleSVD, GetStatisticsMultiscaleSVD, EstimateDim_LargestMultiscaleGap
%

% Copyright (c), 2009
% Mauro Maggioni
% Duke University
%

TimingsEstDim_MSVD = cputime;

%% Handle parameters
if nargin<2,
    cOpts = [];
end;
if ~isfield(cOpts,'MAXDIM'),            cOpts.MAXDIM = size(X,2);       end;
if ~isfield(cOpts,'MAXAMBDIM'),         cOpts.MAXAMBDIM = size(X,2);    end;
if ~isfield(cOpts,'verbose'),           cOpts.verbose = false;          end;
if ~isfield(cOpts,'Ptwise'),            cOpts.Ptwise = false;           end;
if ~isfield(cOpts,'PtIdxs'),            cOpts.PtIdxs = 0;               end;
if ~isfield(cOpts,'NumberOfTrials'),    cOpts.NumberOfTrials = 5;       end;
if ~isfield(cOpts,'NetsOpts'),          cOpts.NetsOpts = [];            end;
if ~isfield(cOpts,'UseSmoothedS'),      cOpts.UseSmoothedS = false;     end;
if ~isfield(cOpts,'Deltas'),            cOpts.Deltas = [];              end;
if ~isfield(cOpts,'EnlargeScales'),     cOpts.EnlargeScales = false;    end;
if ~isfield(cOpts,'KeepV'),             cOpts.KeepV = false;            end;
if ~isfield(cOpts,'kNN'),               cOpts.kNN = [];                 end;
if ~isfield(cOpts,'AllPtsFast'),        cOpts.AllPtsFast = false;       end;
if ~isfield(cOpts,'DownSample'),        cOpts.DownSample = true;        end;
if ~isfield(cOpts,'RandomizedSVDThres'),cOpts.RandomizedSVDThres = inf;  end;

% Multiple trials don't make sense if PtIdxs have been provided
if cOpts.PtIdxs~=0,
    cOpts.NumberOfTrials = 1;
end;

[D,N] = size(X);

% Set interal parameters
AutoScalesMinNPerBin = max([10,round(N/50)]);
MinNetDeltaPts = 0;
MinNetPts = 5;
MinNetDeltaLogRadius=0.99999;
SVDDim = cOpts.MAXDIM;

%% Random projection if requested (not the most efficient algorithm, should use fast random projection
q=[];
if cOpts.MAXAMBDIM < D,
    [q,r] = qr(randn(D,cOpts.MAXAMBDIM));
    q = q';
    q = q(1:cOpts.MAXAMBDIM,:);
    X = q*X;
end;

%% Run the algorithm with cOpts.NumberOfTrials different random multiscale nets
if ~cOpts.Ptwise,
    EstDims           = zeros(cOpts.NumberOfTrials,1);
else
    if cOpts.PtIdxs==0,
        EstDimsPtWise     = zeros(cOpts.NumberOfTrials,size(X,2));
    else
        EstDimsPtWise     = zeros(cOpts.NumberOfTrials,length(cOpts.PtIdxs));
    end;
end;

if ~cOpts.AllPtsFast,
    my_nargout = nargout;       % This fixes what seems to be a Matlab bug, which seem to assign nargout=1 within a parfor loop
    
    %% Go parallel with default configuration
    try
        if matlabpool('size')==0,
            matlabpool
        end;
    catch
    end;
    
    %% Perform multiple trials of dimension estimation
    for k = 1:cOpts.NumberOfTrials,
        Nets = [];
        if isempty(cOpts.NetsOpts),
            % Compute the scales
            if isempty(cOpts.Deltas) && isempty(cOpts.kNN),
                Delta = AutoScales(X,struct('MinNperBin',min([min(size(X,1)),AutoScalesMinNPerBin]),'StatType','min'));
                if cOpts.EnlargeScales,
                    Delta = unique([0.75*Delta,Delta,1.25*Delta]);
                end;
            else
                Delta = cOpts.Deltas;
            end;
            
            % Compute the multiscale svd
            NetsOpts = struct( ...
                'Epsilon',1, ...
                'Delta',Delta, ...
                'NumberOfScales',length(Delta), ...
                'SVDDim',min([size(X,1),size(X,2),cOpts.MAXDIM]), ...
                'ApproxSVD',0, ...
                'RandomizedSVDThres', cOpts.RandomizedSVDThres, ...
                'NormalizedSVD',1, ...
                'DecreasingNets',0, ...
                'DownSample',cOpts.DownSample, ...
                'XIsTransposed', 0, ...
                'Verbose', (cOpts.verbose>1), ...
                'KeepV',cOpts.KeepV, ...
                'MinNetDeltaPts',MinNetDeltaPts, ...
                'MinNetPts',MinNetPts, ...
                'MinNetDeltaLogRadius',MinNetDeltaLogRadius, ...
                'FastNNSearcher','nn_search');
            if ~isempty(cOpts.kNN),
                NetsOpts.kNN = cOpts.kNN;
                NetsOpts.NumberOfScales = length(NetsOpts.kNN);
            end;
        else
            NetsOpts = cOpts.NetsOpts;
        end;
        if cOpts.PtIdxs==0,
            TimingsMultiscaleNet(k) = cputime;
            % Constructs multiscale nets
            Nets = FastRoughSetMultiscaleSVD( X, NetsOpts );
            TimingsMultiscaleNet(k) = cputime - TimingsMultiscaleNet(k);
            NetsOpts.NumberOfScales = length(Nets);
            if isempty(cOpts.kNN),
                NetsOpts.Delta = NetsOpts.Delta(1:length(Nets));
            else
                for j = 1:length(Nets),
                    NetsOpts.Delta(j) = max(Nets(j).Delta);
                end;
            end;
            
            % Computes multiscale statistics
            TimingsMultiscaleStats(k) = cputime;
            if ~cOpts.UseSmoothedS
                MSVD_Stats = PointStatisticsMultiscaleSVD( [], 1:size(X,2), Nets, struct('FindLinearRegion',0,'J',false,'JFr',false,'Beta',false,'Volume',false,'ResVar',false,'UseSmoothedS',false,'Delta',NetsOpts.Delta) );
            else
                MSVD_Stats = GetStatisticsMultiscaleSVD( Nets, struct('NumberOfPoints',size(X,2),'NetsOpts',NetsOpts,'Lambda',100,'Iter',20,'FindLinearRegion',0,'J',false,'JFr',false,'Beta',false,'Volume',false,'ResVar',false) );
            end;
            TimingsMultiscaleStats(k) = cputime-TimingsMultiscaleStats(k);
            % Compute gaps and other stats used for dimensionality estimation
            TimingsEstimateDim(k) = cputime;
            lEstDimStats(k)  = EstimateDim_LargestMultiscaleGap( X, struct('Nets',Nets,'NetsOpts',NetsOpts,'Stats',MSVD_Stats), struct('Verbose',cOpts.verbose,'Ptwise',cOpts.Ptwise) );
            TimingsEstimateDim(k) = cputime - TimingsEstimateDim(k);
            if ~cOpts.Ptwise,
                EstDims(k)   = lEstDimStats(k).DimEst;
            else
                EstDimsPtWise(k,:) = lEstDimStats(k).DimEst;
            end;
        else
            TimingsMultiscaleNet(k) = 0;
            TimingsMultiscaleStats(k) = cputime;
            [lEstDimStats(k),MSVD_Stats(k)] = GetPointSpectrum( X', cOpts.PtIdxs, NetsOpts );
            TimingsMultiscaleStats(k) = cputime-TimingsMultiscaleStats(k);
            EstDimsPtWise(k,:) = lEstDimStats(k).DimEst;             
            TimingsEstimateDim(k) = 0;
        end;        
        
        if my_nargout>1,
            Stats(k).Delta = NetsOpts.Delta;
            Stats(k).NetsOpts = NetsOpts;            
            Stats(k).S = MSVD_Stats.S;
            Stats(k).S_lbl = MSVD_Stats.S_lbl;
            if ~cOpts.PtIdxs,
                Stats(k).Nets = Nets;                
                if cOpts.UseSmoothedS,
                    Stats(k).Tychonov = MSVD_Stats.Tychonov;
                end;
            end;
            if ~isempty(q),
                Stats(k).Q = q;
            end;
        end;
    end;
    
    %% Summarize the estimated dimension from the different trials
    if (~cOpts.Ptwise),
        lGoodIdxs = find(EstDims==round(EstDims));
        if ~isempty(lGoodIdxs),
            EstDim = mean(EstDims(lGoodIdxs));
        else
            EstDim = mean(EstDims);
        end;
        % Summarize EstDimStats
        if my_nargout>1,
            for i = cOpts.NumberOfTrials:-1:1,
                lMinGoodScales(i) = Stats(i).Delta(min(lEstDimStats(i).GoodScales));
                lMaxGoodScales(i) = Stats(i).Delta(max(lEstDimStats(i).GoodScales));
            end;
            EstDimStats.GoodScales = [median(lMinGoodScales),median(lMaxGoodScales)];
            EstDimStats.EstDim     = EstDims;
        end;
    else
        EstDim = median(EstDimsPtWise,1);
        % Summarize EstDimStats
        if my_nargout>1,
            for i = cOpts.NumberOfTrials:-1:1,
                lMinGoodScales(i,:) = Stats(i).Delta(cellfun(@min,lEstDimStats(i).GoodScales));
                lMaxGoodScales(i,:) = Stats(i).Delta(cellfun(@max,lEstDimStats(i).GoodScales));
            end;
            for p = 1:size(lMinGoodScales,2),
                EstDimStats.GoodScales{p} = [median(lMinGoodScales(:,p)),median(lMaxGoodScales(:,p))];
            end;
        end;
    end;
    
    %% If dimension was reduced by random projection, undo the projection to the V's
    if (~isempty(q)) && (isfield(Stats(1).Nets(1),'V') && ~isempty(Stats(1).Nets(1).V)),
        for k = 1:length(Stats),
            for i = 1:length(Stats(k).Nets),
                newV = zeros(size(Stats(k).Nets(i).V,1),D,size(Stats(k).Nets(i).V,3));
                for j = 1:size(Stats(k).Nets(i).V,1),
                    newV(j,:,:) = q*squeeze(Stats(k).Nets(i).V(j,:,:));
                end;
                Stats(k).Nets(i).V = newV;
            end;
        end;
    end;
else
    % Save in appropriate format
    SaveDataPts(X,'EstDim_MSVD_tmp.pts');    
    % Call C program
    [lRet,lOutput]=system('FastSpectralGraph -df EstDim_MSVD_tmp.pts');    
    % Load results
    S = Read3MatrixFromFile('FastSpectralGraph.out');
    
end;

EstDimStats.Timings.EstDim_MSVD       = cputime-TimingsEstDim_MSVD;
EstDimStats.Timings.MultiscaleNet     = mean(TimingsMultiscaleNet);
EstDimStats.Timings.MultiscaleStats   = mean(TimingsMultiscaleStats);
EstDimStats.Timings.EstimateDim       = mean(TimingsEstimateDim);

return;