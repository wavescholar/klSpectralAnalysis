clear all;
EstDimOpts = struct('NumberOfTrials',15,'verbose',1,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',false, 'EnlargeScales',true );
XName = 'CernoTiles'; XNickName = 'C';


figPath = 'K:\AnalysisResults_070413\MediumSet_360\figs\';
mkdir(figPath);

cd('K:\AnalysisResults_070413\MediumSet_360\Tiles_Panel5\V3_0_1_Panel5\');
metadataFile ='73LargestBlobFromEachSlide_Panel5.csv'
versionString = 'V3 0 1 MediumSet 360 Panel5';

fid = fopen(metadataFile);
MetadataFields = { {'basefileName'},  {'nucFile'},  {'marker488'},  {'marker555'},  {'marker647'},  {'Panel'},  {'SlideFileName'},  {'CaseId'},  {'CasePart'},  {'Dx'},  {'Risk'},  {'PatientID'},  {'SurvivalTime'},  {'AgeAtStartCase'}};
DataHeaders =textscan(fid,'%s',1,'delimiter','\n');
DataHeaders = cell2mat(DataHeaders{1});

blobMetadata = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %d %d','delimiter', ',');

numBlobs = size(blobMetadata{1},1);
BlobData = {};
risk = {};

for i=1:numBlobs
    try        
    fileName =[blobMetadata{1,1}{i},'_BlobData.mat'];
    risk{i} = blobMetadata{1,11}{i};
    ldata = load(fileName);
    numBlobs = size(BlobData,2);
    BlobData{numBlobs+1} = ldata;
    catch
        display(['Problem With ',num2str(i)]);
    end
end

%Put All the tiles together
allTileData={};
numAllTiles =0;
riskClass=[];
for i = 1:size(BlobData,2)
    
    tileData = BlobData{i}.tileData;
    [ntd,mtd]=size(tileData);
    
    sr = size(riskClass,2);
    lr = risk{i};
    if strcmp(lr,'high') |  strcmp(lr,'high_recur') | strcmp(lr,'na')
        riskClass(sr+1:sr+ntd)=2;
    end
    if strcmp(lr,'int')
        riskClass(sr+1:sr+ntd)=1;
    end
    if strcmp(lr,'low')
        riskClass(sr+1:sr+ntd)=0;
    end      
        
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

doPCA=0;
if doPCA
    [p c v] = princomp(data_mean);
    X = c(:,1:8);
else
    X =  data_mean;
end

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

for i = 2:1:30    
    E1to3=[Data.G.EigenVecs(:,i),Data.G.EigenVecs(:,i+1),Data.G.EigenVecs(:,i+2)];
    E1to3_High= E1to3(riskClass==2,:);
    E1to3_Int= E1to3(riskClass==1,:);
    E1to3_Low= E1to3(riskClass==0,:);
    scrsz = get(0,'ScreenSize');        
    figure('Visible','off','Position',[1 scrsz(4) scrsz(3) scrsz(4)]); hold on;
    rsH =randsample(size(E1to3_High,1),1200);
    rsI =randsample(size(E1to3_Int,1),1200);
    rsL =randsample(size(E1to3_Low,1),1200);
    scatter3(E1to3_High(rsH,1),E1to3_High(rsH,2),E1to3_High(rsH,3),'r','.');
    scatter3(E1to3_Int(rsI,1),E1to3_Int(rsI,2),E1to3_Int(rsI,3),'y','.');
    scatter3(E1to3_Low(rsL,1),E1to3_Low(rsL,2),E1to3_Low(rsL,3),'b','.');
    title( {['EV ',num2str(i),' ' ,num2str(i+1),' ',num2str(i+2)],versionString})
    legend('High Risk','Int Risk','Low Risk');
    print(gcf,'-dpdf',[figPath,versionString,'EV ',num2str(i),' ' ,num2str(i+1),' ',num2str(i+2),'_0deg.pdf']);
    camroll(45)
    print(gcf,'-dpdf',[figPath,versionString,'EV ',num2str(i),' ' ,num2str(i+1),' ',num2str(i+2),'_45deg.pdf']);
    camroll(135)
    print(gcf,'-dpdf',[figPath,versionString,'EV ',num2str(i),' ' ,num2str(i+1),' ',num2str(i+2),'_135deg.pdf']);
    close gcf; delete gcf;
end


for i = 2:1:30%i = 2:3:30    
    E1to3=[Data.G.EigenVecs(:,i),Data.G.EigenVecs(:,i+1)];
    E1to3_High= E1to3(riskClass==2,:);
    E1to3_Int= E1to3(riskClass==1,:);
    E1to3_Low= E1to3(riskClass==0,:);
    scrsz = get(0,'ScreenSize');        
    figure('Visible','off','Position',[1 scrsz(4) scrsz(3) scrsz(4)]); hold on;
    rsH =randsample(size(E1to3_High,1),1200);
    rsI =randsample(size(E1to3_Int,1),1200);
    rsL =randsample(size(E1to3_Low,1),1200);
    scatter(E1to3_High(rsH,1),E1to3_High(rsH,2),'r','.');
    scatter(E1to3_Int(rsI,1),E1to3_Int(rsI,2),'y','.');
    scatter(E1to3_Low(rsL,1),E1to3_Low(rsL,2),'b','.');
    title( {['EV ',num2str(i),' ' ,num2str(i+1)],versionString})
    legend('High Risk','Int Risk','Low Risk');
    print(gcf,'-dpdf',[figPath,versionString,'2DEV ',num2str(i),' ' ,num2str(i+1),'_.pdf']);
    close gcf; delete gcf;
end



