function Data = Analyse_MSVD_Data( cDataName, cOpts )

%
% function Data = Analyse_MSVD_Data( cDataName, cOpts )
%
% Function for diffusion and multiscale svd analysis of data. Global variables:
%   X   : N by D matrix of N points in R^D
%
% IN:
%   cDataName   :   name of MD data. Used to load data files, see below
%   [cOpts]     :   structure with the following fields:
%                   [RestrictIdxs]  : indices to which to restrict the analysis.
%                   [GraphDiffOpts] : structure of options to be passed to GraphDiffusion. 
%                                     If it's a string, loads from filename. Default: <cDataName>\<cDataName>_GraphDiffOpts.mat,
%                                     if not found, defaults to a set of standard options.
%                   [NetsOpts]      : structure of options to be passed to FastRoughSetMultiscaleSVD. 
%                                     If it's a string, loads from filename. Default: <cDataName>\<cDataName>_NetsOpts.mat,
%                                     if not found, defaults to a set of standard options.
%                   [SaveFileName]  : filename to save the results to. If empty, does not save results. 
%                                     Default: <cDataName>\<cDataName>_Data.mat
%                   [MSVD]          : perform multiscale SVD analysis. Default: true.
%               
%
% OUT:
%   Data        :   structure with the following fields:
%                   G       : structure as returned by GraphDiffusion
%                   Nets    : multiscale nets as returned by FastRoughSetMultiscaleSVD
%                   NetsOpts: structure of options used to compute the multiscale svd.
%                   NetsOpts.DistInfo : matrix of distances used by FastRoughSetMultiscaleSVD, 
%                               returned only if smaller than that in <cDataName>_DistInfo.
%                   Stats   : structure of multiscale statistics as returned by GetStatisticsMultiscaleSVD
%                   Coords  : M by N matrix of M configurations diffusion-mapped into R^N
%                   DataViewerOpts : structure of options useful for calling DataViewer
%

% USES:
%   Files       :   <cDataName>\<cDataName>_X, <cDataName>\<cDataName>_DistInfo
%   Functions   :   NormalizeMolecule, GraphDiffusion, RestrictDistInfoToSubSet, 
%                   FastRoughSetMultiscaleSVD, GetStatisticsMultiscaleSVD

%
%
% (c) Copyright 
% Mauro Maggioni
% Duke University, 2008
%

global X DistInfo
global Data                         % This is for debugging purposes/not recompute things that already exist

% Load distance data if needed
fprintf('Loading data...');
if isempty(X),
    % This is the file with all configurations
    load(sprintf('%s\\%s_X',cDataName,cDataName));
end;
% This is the file with all distances
if isempty(DistInfo),
    lFileName = sprintf('%s\\%s_DistInfo',cDataName,cDataName);
    if exist(lFileName),
        load(lFileName);
    else
        DistInfo = [];
    end;
end;
fprintf('done.\n');

if nargin<2,
    cOpts = [];
end;
if ~isfield(cOpts,'RestrictIdxs'),  cOpts.RestrictIdxs = [];    end;
if ~isfield(cOpts,'GraphDiffOpts'), cOpts.GraphDiffOpts= [];    end;
if ~isfield(cOpts,'NetsOpts'),      cOpts.NetsOpts= [];         end;
if ~isfield(cOpts,'SaveFileName'),  cOpts.SaveFileName = sprintf('%s\\%s_Data.mat',cDataName,cDataName); end;
if ~isfield(cOpts,'MSVD'),          cOpts.MSVD = true;          end;

if ~isempty(DistInfo),
    if ~isfield(DistInfo,'isinidxs'),
        fprintf('Sorting indices and distances for faster nn searches...');
        AddSortedIdxsToDistInfo;
        fprintf('done.\n');
    end;
end;


%%
%
% Handling options
%
%
if ~isempty(cOpts.RestrictIdxs),
    fprintf('Restricting data and distances to subset...');
    X = X(cOpts.RestrictIdxs,:);
    if ~isempty(DistInfo),
        DistInfo = RestrictDistInfoToSubSet( DistInfo, cOpts.RestrictIdxs );
    end;
    fprintf('done.\n');
end;

GraphDiffOpts = [];
if ~isempty(cOpts.GraphDiffOpts),
    if isstr(cOpts.GraphDiffOpts),
        if exist(cOpts.GraphDiffOpts),
            load(cOpts.GraphDiffOpts); end;
    elseif isstruct(cOpts.GraphDiffOpts),
        GraphDiffOpts = cOpts.GraphDiffOpts; end;
end;
if isempty(GraphDiffOpts),
    if exist(sprintf('%s\\%s_GraphDiffOpts.mat',cDataName,cDataName)),
        load(sprintf('%s\\%s_GraphDiffOpts.mat',cDataName,cDataName)); end;
    if isempty('GraphDiffOpts'),
        % Set default options for constructing the graph
        GraphDiffOpts.Normalization = 'sFokkerPlanck';
        GraphDiffOpts.Epsilon = 1;
        GraphDiffOpts.kNN = 100;
        GraphDiffOpts.kNNAutotune = 20;
        GraphDiffOpts.kEigenVecs = 35;
        GraphDiffOpts.Symmetrization = 'W+Wt';
        GraphDiffOpts.DontReturnDistInfo = 1;
    end;
end;
GraphDiffOpts.DistInfo = DistInfo;                      %GraphDiffOpts.Metric = @MolecularMetric;

NetsOpts = [];
if ~isempty(cOpts.NetsOpts),
    if isstr(cOpts.NetsOpts),
        if exist(cOpts.NetsOpts),
            load(cOpts.NetsOpts); end;
    elseif isstruct(cOpts.NetsOpts),
        NetsOpts = cOpts.NetsOpts; end;
end;
if isempty(NetsOpts),
    if exist(sprintf('%s\\%s_NetsOpts.mat',cDataName,cDataName)),
        load(sprintf('%s\\%s_NetsOpts.mat',cDataName,cDataName)); end;
    if isempty(NetsOpts),
        NetsOpts = struct('Epsilon',1,'Delta',[2:0.25:10],'NumberOfScales',length([2:0.25:10]), ...
            'SVDDim',10,'ApproxSVD',0,'NormalizedSVD',1,'DecreasingNets',0, ...
            'DownSample',1,'KeepV',0,'MinNetPts',30);
    end;
end;
if ~isfield(Data,'NetsOpts'),
    Data.NetsOpts = NetsOpts; clear NetsOpts;
    Data.NetsOpts.DistInfo = DistInfo;
end;



%%
%
% Some preprocessing
%
%

% Normalize the molecules by center of mass - TODO: need mass information, H information etc...
%Data.X_norm = NormalizeMolecule( X );
lNumberOfPoints = size(X,1);


%%
%
% DIFFUSION MAP, with self-tuning
%
%

if (~isfield(Data,'G')) || (isempty(Data.G)),
    fprintf('\n\nConstructing graph and diffusion map...\n');
    Data.G = GraphDiffusion(X, 0, GraphDiffOpts);                  % These are for debugging purposes
    fprintf('done.\n');
end;

Data.Coords = Data.G.EigenVecs;





%%
%
% MULTISCALE SVD
%
%
if cOpts.MSVD,
    if (~isfield(Data,'Nets')) || (isempty(Data.Nets)),                 % These are for debugging purposes
        fprintf('\n\n Performing Multiscale SVD analysis...\n');
        Data.Nets = FastRoughSetMultiscaleSVD( X, Data.NetsOpts );
        Data.NetsOpts.NumberOfScales = length(Data.Nets);
        Data.NetsOpts.Delta = Data.NetsOpts.Delta(1:length(Data.Nets));
        fprintf('done.\n');
    end;
    
    % Computes statistics
    if (~isfield(Data,'Stats')),
        fprintf('Computing Multiscale Statistics...\n');
        Data.Stats = GetStatisticsMultiscaleSVD( Data.Nets, struct('NumberOfPoints',lNumberOfPoints,'NetsOpts',Data.NetsOpts,'Lambda',40,'Iter',20,'J',true,'JFr',true,'Volume',true,'Beta',true,'ResVar',true) );
        fprintf('done.\n');
    end;
end;

%
% Save file
%
if ~isempty(cOpts.SaveFileName),
    if exist(cOpts.SaveFileName),
        reply = input(sprintf('Output file %s already exists. Overwrite? Y/N [Y]: ',cOpts.SaveFileName), 's');
        if ~(strcmpi(reply,'y') || strcmpi(reply,'yes')),
            return;
        end;
    end;
    fprintf('\n Saving...');
    if isempty(cOpts.RestrictIdxs),
        Data.NetsOpts.DistInfo=[]; end;
    eval(['save ',cOpts.SaveFileName,' -struct Data']);
end;
fprintf('done.\n');

return;

