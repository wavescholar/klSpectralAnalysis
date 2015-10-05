%% Clear variables
clear all;close all;clc; clear global

%% Generate data set
fprintf('\n Generating data...');
[X0,~,~,Labels]=GenerateDataSets('Planes',struct('NumberOfPoints',1000,'EmbedDim',3,'PlanesOpts',struct('Number',2,'Dim',[1,2],'Origins',zeros(2,3),'NumberOfPoints',[100,1900],'Dilations',[1,5])));
X=X0(find(Labels==1),:);
Y=X0(find(Labels==2),:);
fprintf('done.');

%% Generate integral operator
fprintf('\n Generating operator...');
K=GenerateIntegralOperatorFromSets(X,Y)';
fprintf('done.');

%% Compute geometric wavelets
fprintf('\n Computing GWT and transforming data...');
GWTopts = struct('pruning',3,'precision',1e-4);
GWT = GMRA(K);
Data = FGWT_trainingData(GWT, K);
%Data.CelWavCoeffs(:,2:end) = threshold_coefficients(Data.CelWavCoeffs(:,2:end), struct('coeffs_threshold',GWTopts.precision,'shrinkage','hard));
[GWT, Data] = simplify_the_GWT_tree(GWT, Data);
[Data.Projections, Data.TangentialCorrections] = IGWT_trainingData(GWT, Data.CelWavCoeffs);
fprintf('done.');

%% Figures
fprintf('\n Generating figures...');
figure;treeplot(GWT.cp);
for j = 1:max(GWT.Scales),
    %figure;scatter3(Y(:,1),Y(:,2),Y(:,3),20,get_partition_at_scale(j,GWT.Scales,GWT.PointsInNet),'filled');colormap;
    figure;scatter3(Y(:,1),Y(:,2),Y(:,3),20,get_partition_at_scale(GWT,j),'filled');colormap;
    title(sprintf('Partition at scale %d.',j));
end;

fprintf('done.');

%% Scaling functions and wavelet display and compression
[lSortedLinePts,lSortedLinePtsIdxs]=sort(X(:,1));
ScalFuns=cat(2,GWT.ScalFuns{:});
WavFuns=cat(2,GWT.WavBases{:});
for k = 1:size(ScalFuns,2),if sum(ScalFuns(:,k).*ScalFuns(:,1))<0,ScalFuns(:,k)=-ScalFuns(:,k);end;end;
for k = 1:size(WavFuns,2),if sum(WavFuns(:,k).*WavFuns(:,1))<0,WavFuns(:,k)=-WavFuns(:,k);end;end;
figure;plot(ScalFuns(lSortedLinePtsIdxs,:));title('All scaling functions');
figure;plot(WavFuns(lSortedLinePtsIdxs,:));title('All wavelet functions');
lMeanScalFuns = mean(ScalFuns,2);
[U,S,V] = svd(ScalFuns-repmat(lMeanScalFuns,[1,size(ScalFuns,2)]));S=diag(S);
figure;plot(log10(S));title('SVD of scaling functions');
lMeanWavFuns = mean(WavFuns,2);
[U,S,V] = svd(WavFuns-repmat(lMeanWavFuns,[1,size(WavFuns,2)]));S=diag(S);
figure;plot(log10(S));title('SVD of wavelet functions');
fprintf('\n');


%% Work on the transpose
Xo=X;
Yo=Y;
X=X';
Kt=K';
GWTt = GMRA(Kt);
[GWTt, Datat] = FGWT_trainingData(GWTt, Kt);
%Data.CelWavCoeffs(:,2:end) = threshold_coefficients(Data.CelWavCoeffs(:,2:end), struct('coeffs_threshold',GWTopts.precision,'shrinkage','hard));
[GWTt, Datat] = simplify_the_GWT_tree(GWTt, Datat);
[Datat.Projections, Datat.TangentialCorrections] = IGWT_trainingData(GWTt, Datat.CelWavCoeffs);

GWTo=GWT;
GWT=GWTt;

%% Figures
fprintf('\n Generating figures...');
figure;treeplot(GWT.cp);
for j = 1:max(GWT.Scales),
    figure;scatter3(X(:,1),X(:,2),X(:,3),20,get_partition_at_scale(j,GWT.Scales,GWT.PointsInNet),'filled');colormap;
    title(sprintf('Partition at scale %d.',j));
end;

fprintf('done.');

%% Scaling functions and wavelet display and compression
[lSortedLinePts,lSortedLinePtsIdxs]=sort(X(:,1));
ScalFuns=cat(2,GWT.ScalFuns{:});
WavFuns=cat(2,GWT.WavBases{:});
for k = 1:size(ScalFuns,2),if sum(ScalFuns(:,k).*ScalFuns(:,1))<0,ScalFuns(:,k)=-ScalFuns(:,k);end;end;
for k = 1:size(WavFuns,2),if sum(WavFuns(:,k).*WavFuns(:,1))<0,WavFuns(:,k)=-WavFuns(:,k);end;end;
figure;plot(ScalFuns(lSortedLinePtsIdxs,:));title('All scaling functions');
figure;plot(WavFuns(lSortedLinePtsIdxs,:));title('All wavelet functions');
lMeanScalFuns = mean(ScalFuns,2);
[U,S,V] = svd(ScalFuns-repmat(lMeanScalFuns,[1,size(ScalFuns,2)]));S=diag(S);
figure;plot(log10(S));title('SVD of scaling functions');
lMeanWavFuns = mean(WavFuns,2);
[U,S,V] = svd(WavFuns-repmat(lMeanWavFuns,[1,size(WavFuns,2)]));S=diag(S);
figure;plot(log10(S));title('SVD of wavelet functions');
fprintf('\n');

GWT=GWTo;

return;



