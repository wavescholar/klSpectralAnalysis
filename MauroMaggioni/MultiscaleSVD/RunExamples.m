%
% This file runs several examples and creates some of the figures in the paper
% "Multiscale Geometric Methods for Data Sets I: Intrinsic dimension and L^2 curvatures", A. Little, M. Maggioni, L. Rosasco, 2011
%
% Execute 'clear all' before running the script, otherwise it will re-run with any existing options without prompting the user.
%

pExampleNames   = { '9-d sphere with little noise', ...
    'Figure 1a,b in [LMR] paper (9-d sphere with large noise)', ...
    'Figure 1c in [LMR] paper (9-d sphere with large noise, single point)', ...
    'Figure 8a,b in [LMR] paper (sphere and segment)', ...
    'Figure 8c,d in [LMR] paper (spiral and plane)', ...
    '47-dimensional sphere, 8000 points, large noise', ...
    'Meyerstaircase', ...
    '6-D Cube with low sampling rate', ...
    'S-Manifold', ...
    '10-D cube with high noise', ...
    '9-D sphere with high noise', ...
    'Two lines and a plane', ...
    'Isomap Faces',...
    'CBCLFaces1',...
    'CBCLFaces2',...
    'HandVideo',...
    'FaceVideo',...
    'ScienceNews', ...
    'Patches from image of Lena'};

fprintf('\n\n Select example to run:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

while true,
    if (~exist('pExampleIdx') || isempty(pExampleIdx) || pExampleIdx==0),
        try
            pExampleIdx = input('');
            pExampleIdx = str2num(pExampleIdx);
        catch
        end;
    end;
    if (pExampleIdx>=1) && (pExampleIdx<=length(pExampleNames)),
        break;
    else
        fprintf('\n %d is not a valid Example. Please select a valid Example above.',pExampleIdx);
        pExampleIdx=0;
    end;
end;

%% Set parameters for dimensionality estimation
EstDimOpts = struct('NumberOfTrials',1,'verbose',0,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',true, 'EnlargeScales',true );   % Usually set NumberOfTrials=10

%% Set parameters for data generation and generates the data set
switch pExampleIdx
    case 1
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',2000,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01/sqrt(100));
    case 2
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
    case 3
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        EstDimOpts.Ptwise = true;
        EstDimOpts.PtIdxs = 1;
    case 4
        XName = 'SphereAndLine'; XNickName = 'Sphere and Line';
        XOpts = struct('Dim','2.5','NoiseType','Gaussian','NoiseParam',0.00);
        EstDimOpts.Ptwise = true;
    case 5
        XName = 'SpiralAndPlane'; XNickName = 'Spiral and Plane';
        XOpts = struct('NumberOfPoints',1100,'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels=[ones(300,1);2*ones(800,1)];
    case 6
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',4000,'Dim',47,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01);
    case 7
        XName = 'MeyerStaircase'; XNickName = 'Z';
        XOpts = struct('NumberOfPoints',500,'Dim',1000,'MeyerStepWidth',20,'EmbedDim',1000,'NoiseType','Gaussian','NoiseParam',0.05/sqrt(1000));
    case 8
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',250,'Dim',6,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01/sqrt(100));
    case 9
        XName = 'S-Manifold'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01);
    case 10
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',1000,'Dim',10,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
    case 11
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        EstDimOpts.EnlargeScales = true;
    case 12
        XName = 'TwoLinesAndAPlane'; XNickName = 'Two Lines and a Plane';
        XOpts = struct('NumberOfPoints',800,'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels= [ones(400,1);2*ones(400,1)];
    case 13
        XName = 'Isomapfaces'; XNickName = 'Isomap Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 14
        XName = 'CBCLFaces1'; XNickName = 'CBCLFaces1 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 15
        XName = 'CBCLFaces2'; XNickName = 'CBCLFaces2 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 16
        XName = 'HandVideo'; XNickName = 'FaceVideo Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 17
        XName = 'FaceVideo'; XNickName = 'FaceVideo Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 18
        XName = 'ScienceNews'; XNickName = 'ScienceNews Faces';
        XOpts = struct(); %struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 19
        XName = 'ImagePatches'; XNickName = 'Lena Patches';
        XOpts = struct('ImageFileName','Lena.jpg','PatchSize',16,'DownsampleImage',8);
        EstDimOpts.Ptwise = true;
end;

%% Generate the data sets
fprintf('\nGenerating %s data...', XName);
[X,GraphDiffOpts,NetsOpts] = GenerateDataSets( XName, XOpts);
X = 1*X;
fprintf('done.');

%% Do Multiscale SVD analysis and estimate intrinsic dimension.
EstDimOpts.RandomizedSVDThres = min([size(X,1),size(X,2),100]);      %EstDimOpts.kNN = round(size(X,1)/100:max([size(X,1)/100,1]):size(X,1)/2);
fprintf('\nEstimating %s dimensionality...', XName );
[EstDim,EstDimStats,Stats] = EstDim_MSVD( X, EstDimOpts );
fprintf('done.');
if ~EstDimOpts.Ptwise
    fprintf('\n\nFinal dimension estimate: %d',EstDim);
    if EstDim~=round(EstDim),
        fprintf('\n Dimension estimates during various randomized runs:'); fprintf('%d,',EstDimStats.EstDim);
    end;
    fprintf('\n\n');
else
    figure;plot(sort(EstDim));axis tight; title('Estimate ptwise dimensionality');
end;

%% Generate plot of multiscale singular values and their standard deviation
if ~EstDimOpts.Ptwise,
    for k = 1:length(EstDimStats),
        if length(EstDimStats(k).GoodScales)>3,
            break;
        end;
    end;
    figure;
    for j = 1:size(Stats(1).S,3),
        PlotWithStd(Stats(1).Delta,(squeeze(Stats(1).S(:,:,j))));hold on;
    end; axis tight;
    YLims = [min(min(mean(Stats(1).S,2))),max(max(mean(Stats(1).S,2)))];set(gca,'YLim',[YLims]);
    % MM: <patch> has bugs that prevent it from producing correct displays if the scale of the numbers involved is very small.
    if ~isempty(EstDimStats(k).GoodScales),
        p=patch([min(EstDimStats.GoodScales),min(EstDimStats.GoodScales),max(EstDimStats.GoodScales),max(EstDimStats.GoodScales)],[min(YLims),max(YLims),max(YLims),min(YLims)],'g');set(p,'FaceAlpha',0.1);
    end;
    title(sprintf('%s^{%d}(%d,%d,%.2f): multiscale average singular values and their standard deviations',XNickName, XOpts.Dim,XOpts.NumberOfPoints,XOpts.EmbedDim,XOpts.NoiseParam));
    xlabel('r'); ylabel('E_z[\sigma_i(z,r)]');
    ylim = get(gca,'YLim'); set(gca,'YLim',[0,max(ylim)]);                  %print(gcf,'-depsc2','CodeForPaper/Figures/Figure1_MSVD_plot.eps');
    
    if isfield(Stats(1),'Tychonov'),
        figure;
        for j = 1:size(Stats(1).S,3),
            PlotWithStd(Stats(1).Delta,(squeeze(Stats(1).Tychonov.S(:,:,j))));hold on;
        end; axis tight;
        YLims = [min(min(mean(Stats(1).Tychonov.S,2))),max(max(mean(Stats(1).Tychonov.S,2)))];set(gca,'YLim',[YLims]);
        if ~isempty(EstDimStats(k).GoodScales),
            p=patch([min(EstDimStats.GoodScales),min(EstDimStats.GoodScales),max(EstDimStats.GoodScales),max(EstDimStats.GoodScales)],[min(YLims),max(YLims),max(YLims),min(YLims)],'g');set(p,'FaceAlpha',0.1);
        end;
        title(sprintf('%s^{%d}(%d,%d,%.2f): multiscale average smoothed singular values and their standard deviations',XNickName, XOpts.Dim,XOpts.NumberOfPoints,XOpts.EmbedDim,XOpts.NoiseParam));
        xlabel('r'); ylabel('E_z[\sigma_i(z,r)]');
        ylim = get(gca,'YLim'); set(gca,'YLim',[0,max(ylim)]);
    end;
elseif length(EstDim)==1,
    figure;
    plot(Stats(1).Delta,squeeze(Stats(1).S),'b');
    YLims = [min(min(mean(Stats(1).S,2))),max(max(mean(Stats(1).S,2)))];set(gca,'YLim',[YLims]);
    if ~isempty(EstDimStats(1).GoodScales),
        p=patch([min(EstDimStats.GoodScales{1}),min(EstDimStats.GoodScales{1}),max(EstDimStats.GoodScales{1}),max(EstDimStats.GoodScales{1})],[min(YLims),max(YLims),max(YLims),min(YLims)],'g');set(p,'FaceAlpha',0.1);
    end;
    title(sprintf('%s^{%d}(%d,%d,%.2f): multiscale singular values at z',XNickName, XOpts.Dim,XOpts.NumberOfPoints,XOpts.EmbedDim,XOpts.NoiseParam));
    xlabel('r'); ylabel('\sigma_i(z,r)');
    ylim = get(gca,'YLim'); set(gca,'YLim',[0,max(ylim)]);
end;

% If pointwise estimation, plot dimension estimation as a function on the data
if EstDimOpts.Ptwise && length(EstDim)==size(X,1),
    figure;scatter3(X(1,:),X(2,:),X(3,:),20,EstDim,'filled');colorbar;title('Pointwise dimension estimate');
    figure;scatter3(X(1,:),X(2,:),X(3,:),20,round(EstDim),'filled');colorbar;title('Pointwise integral dimension estimate');
end;

% Plot maximal good scale as a function of the point
if ~EstDimOpts.Ptwise,
    figure;scatter3(X(1,:),X(2,:),X(3,:),20,repmat(max(EstDimStats.GoodScales),[size(X,2),1]),'filled');colorbar;title('Maximal good local scale');
elseif length(EstDimStats.GoodScales)==size(X,1),
    figure;scatter3(X(1,:),X(2,:),X(3,:),20,cellfun(@max,EstDimStats.GoodScales),'filled');colorbar;title('Maximal good local scale');
end;

% Generate other plots
MGaps = -diff(Stats(1).S,1,3);
MGaps_mean = squeeze(mean(MGaps,2));
figure;imagesc(flipdim(MGaps_mean(:,1:min([4*ceil(EstDim),size(MGaps_mean,2)])),1));colorbar;set(gca,'YTickLabel',fliplr(Stats(1).Delta(max(1,floor(get(gca,'YTick'))))));xlabel('Sing. Val.');ylabel('Scale');title('Multiscale Gaps');grid on;colormap(gray);

if isfield(Stats(1),'Tychonov'),
    MGapsT = -diff(Stats(1).Tychonov.S,1,3);
    MGaps_meanT = squeeze(mean(MGapsT,2));
    % Surface plot
    z=double(flipdim(MGaps_meanT(:,2:min([ceil(EstDim)+5,size(MGaps_meanT,2)])),1));
    x=(1:size(z,2))+1;y=Stats(1).Delta(1:size(z,1));
    [xi,yi]=meshgrid(linspace(min(x),max(x),50),linspace(min(y),max(y),50));zi=interp2(x,y,z,xi,yi,'linear');
    figure;surf(xi,yi,zi,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');axis tight;hold on;
    %[xi,yi]=meshgrid(x,y);
    p=plot3(xi,yi,zi,'k.','MarkerSize',4);xlabel('r');ylabel('E_z[\sigma_i(z,r)]');alpha(0.7);grid off;
    title('Multiscale Gaps as a function of scale and gap number');
end;

try
    Data.Delta = Stats(1).Delta;
    Data.NetsOpts = Stats(1).NetsOpts;
    Data.Stats = Stats(1);
    Data.Nets = Stats(1).Nets;
catch
end;

if pExampleIdx==5,
    lEstDim = round(EstDim)';                % Round off dimensionality
    lEstDim(lEstDim>2)=2;                    % Since we know there are two clusters (that's an input to TPMM), and 1 and 2 are the most frequent estimated dimensions, we arbitrarily assign any point of estimated dimension >2 to the cluster of dimension 2. In fact, we know more since the algorithm correctly identifies the intersection points.
    fprintf('\n Classification rate: %.2f\%]n',sum(round(lEstDim)==Labels)/length(Labels));
end;

fprintf('\n');

%% Now create things for the UI
global X Data

Data.EstDim = EstDim;
Data.EstDimStats = EstDimStats;
Data.XCoords = X';
Data.Coords = X';

fprintf('\n\nConstructing graph and diffusion map...\n');
Data.G = GraphDiffusion(X, 0, GraphDiffOpts);                  % These are for debugging purposes
fprintf('done.\n');


if (exist('Labels')) && (~isempty(Labels)),
    Data.Coords = [Data.Coords,Labels];
    Data.G.EigenVecs=[Data.G.EigenVecs,Labels];
end;

fprintf('\n');
X=X';
global X Data
GUI_Analyse_MD



return;