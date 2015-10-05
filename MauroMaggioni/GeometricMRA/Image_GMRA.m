function [gMRA,gMRA_Data]=Image_GMRA( Image, opts )

% function Image_GMRA( Image, opts )
%
% IN:
%   Image   : an image, as either a matrix or a string, indicating the image name (compatible with being loaded by imread)
%   [opts]  : structure of options:
%               [PatchSize]   : size of patches. Default: 16.
%               [Downsample]  : downsampling factor for image. Default: 1 (no downsampling).
%       
%
% OUT:
%   gMRA        : geometric wavelets structure
%   gMRA_Data   : GWT of the data
% 
% EXAMPLE:
%   [gMRA,gMRA_Data]=Image_GMRA( 'Lena.jpg', struct('PatchSize',16,'Downsample',4) );
%

% @Copyright
% Mauro Maggioni
% Duke University, 2010
%

%% Parameters
if nargin<2,                        opts = [];              end;
if ~isfield(opts,'PatchSize'),      opts.PatchSize = 16;    end;
if ~isfield(opts,'Downsample'),     opts.Downsample = 1;    end;

%% Load image
if ischar(Image),
    I = imread(Image);    
else
    I = Image;
end;

I = double(I);
I=I/max(max(I));

%% Downsample image
if opts.Downsample>1,
    I = I(1:opts.Downsample:size(I,1),1:opts.Downsample:size(I,2));
end;

%% Construct data set
[gMRA.X,XCoords] = filtergraph(I,'ptch',opts.PatchSize); 


%% set GWT parameters
GWTopts = struct();
GWTopts.ManifoldDimension = 0; % if 0, then determine locally adaptive dimensions using the following fields:
GWTopts.errorType = 'relative'; % or absolute
GWTopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
GWTopts.precision  = 1e-3; % only for leaf nodes
% The following thresholds are used in the code construct_GMRA.m
GWTopts.threshold1 = 1e-1; % threshold of singular values for determining the rank of each ( I - \Phi_{j,k} * \Phi_{j,k} ) * Phi_{j+1,k'}
GWTopts.threshold2 = 5e-2; % threshold for determining the rank of intersection of ( I - \Phi_{j,k} * \Phi_{j,k} ) * Phi_{j+1,k'}
% The following parameter .pruning determines which version of geometric wavelets to use
GWTopts.pruning = 1;
% whether to use best approximations
GWTopts.addTangentialCorrections = true;
% whether to sparsify the scaling functions and wavelet bases
GWTopts.sparsifying = false;
% whether to split the wavelet bases into a common intersection and
% children-specific parts
GWTopts.splitting = false;
% METIS parameters
GWTopts.knn = 50;
GWTopts.knnAutotune = 30;
GWTopts.smallestMetisNet = 20;


%% Construct Geometric Wavelets
fprintf('\n Constructing Geometric Wavelets...');
gMRA = GMRA(gMRA.X,GWTopts);
gMRA.XCoords = XCoords;
fprintf('done.');

%% Computing all wavelet coefficients 
fprintf('\n Computing GWT of original data...');
[gMRA,gMRA_Data] = GWT_trainingData(gMRA, gMRA.X);
fprintf('done.');

fprintf('\n');

%% Display results
% Display the GWT coefficients
GWT_DisplayCoeffs( gMRA, gMRA_Data );

% Plot GWTapproximation error
GWT_DisplayApproxErr( gMRA, gMRA_Data );

% Plot some elements in the dictionary
figure;treeplotimages(gMRA,7);

figure;imagesc(reshape(squeeze(gMRA.X(:,10)),size(I)));colormap(gray);title('Original downsampled image');
figure;for k = 1:size(gMRA_Data.Projections,3);imagesc(reshape(squeeze(gMRA_Data.Projections(:,128,k)),size(I)));colormap(gray);title(sprintf('Reconstruction at scale %d',k));pause;end;

gMRA.Image = I;

return;