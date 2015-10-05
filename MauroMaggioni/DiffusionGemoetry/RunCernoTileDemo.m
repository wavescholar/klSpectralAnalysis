EstDimOpts = struct('NumberOfTrials',15,'verbose',1,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',false, 'EnlargeScales',true );
XName = 'CernoTiles'; XNickName = 'C';

cd('d:\klSpectralAnalysis\MauroMaggioni\DiffusionGemoetry\TileData');
files = dir

%Load The Blob Data
BlobData = {};
for i=3:size(files,1)
    fileName = files(i).name;
    ldata = load(fileName);
    numBlobs = size(BlobData,2);
    BlobData{numBlobs+1} = ldata;
end

%Put All the tiles together
allTileData={};
numAllTiles =0;
for i = 1:size(BlobData,2)
    
    tileData = BlobData{i}.caseData;
    [ntd,mtd]=size(tileData);
    numAllTiles = numAllTiles+ntd
    [natd,matd]=size(allTileData);
    
    for j=1:ntd
        allTileData{matd+j}=tileData{j};
    end
end
allTileData =  allTileData'

featureIncludeIndex=[1:117];
 
ntiles = size(allTileData,1); 
nfeats = length(featureIncludeIndex);

data_mean = zeros(ntiles,nfeats);
data_std = zeros(ntiles,nfeats);
data_quantile = zeros(ntiles,nfeats);

for tileCount = 1:ntiles
    tile = allTileData{tileCount};
    for f = 1:size(featureIncludeIndex,2)
        feattile = tile{featureIncludeIndex(f),1};
        if ~isempty(feattile)
            data_mean(tileCount,f) = mean(feattile); 
            data_std(tileCount,f) = std(feattile);
            data_quantile(tileCount,f) = quantile(feattile,.95);
        end
    end
end

[numRows numColumns] = size(data_mean)
%Get rid of any NaN's
for rowI = 1:numRows
    for colI = 1:numColumns
        if isnan(data_mean(rowI,colI))
            disp([sprintf('%f',rowI), ' '  ,sprintf('%f',colI)]);
            data_mean(rowI,colI) = 0;
        end
    end
end

if doPCA
    [p c v] = princomp(data_mean);
    data = c(:,1:8);
else
    data = data_quantile;
end

X = data_mean;

%'EmbedDim',size(X,2) - THIS CAN NOT BE RIGHT

XOpts = struct('NumberOfPoints',size(X,1),'Dim',size(X,2),'EmbedDim',size(X,2),'NoiseType','Gaussian','NoiseParam',0.1);

GraphDiffOpts = struct( ...
    'Normalization','smarkov', ...
    'Epsilon',1, ...
    'kNN', 100, ...
    'kNNAutotune', 20, ...
    'kEigenVecs', 35, ...
    'Symmetrization', 'W+Wt', ...
    'DontReturnDistInfo', 1, ...
    'XIsTransposed', 1);
GraphDiffOpts.DontReturnDistInfo=0;
Data.G = GraphDiffusion(X, 0, GraphDiffOpts);                  % These are for debugging purposes

figure;
subplot(1,2,1);plot3(Data.G.EigenVecs(:,2),Data.G.EigenVecs(:,3),Data.G.EigenVecs(:,4),'.');title('Diffusion embedding (2,3,4)');
subplot(1,2,2);plot3(Data.G.EigenVecs(:,5),Data.G.EigenVecs(:,6),Data.G.EigenVecs(:,7),'.');title('Diffusion embedding (5,6,7)');

figure;hold on;
plot(Data.G.EigenVecs(:,1));
plot(Data.G.EigenVecs(:,2));
plot(Data.G.EigenVecs(:,3));
plot(Data.G.EigenVecs(:,4));
plot(Data.G.EigenVecs(:,5));
plot(Data.G.EigenVecs(:,6));
plot(Data.G.EigenVecs(:,7));
