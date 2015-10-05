clear all;close all;clc; clear global

%% Set parameters
pValN = 3;
pMaxDimForNNs = inf;

DataSet(1).name = 'SwissRoll';
DataSet(1).opts = struct('NumberOfPoints',2000,'EmbedDim',10);
DataSet(2).name = 'BMark_MNIST';
DataSet(2).opts = struct('NumberOfPoints',2000,'MnistOpts',struct('Sampling', 'RandN', 'QueryDigits',7, 'ReturnForm', 'vector'));

% Set parameters for GWT
GWTopts{1} = struct('GWTversion',0);
GWTopts{1}.ManifoldDimension = 2;
GWTopts{1}.threshold1 = 1e-3;
GWTopts{1}.threshold2 = .1;
GWTopts{1}.addTangentialCorrections = false;
GWTopts{1}.sparsifying = false;
GWTopts{1}.splitting = false;
GWTopts{1}.knn = 30;
GWTopts{1}.knnAutotune = 20;
GWTopts{1}.smallestMetisNet = 10;
GWTopts{1}.verbose = 1;
GWTopts{1}.shrinkage = 'hard';
GWTopts{1}.avoidLeafnodePhi = false;
GWTopts{1}.mergePsiCapIntoPhi  = true;
GWTopts{1}.coeffs_threshold = 0;
GWTopts{1}.errorType = 'relative';
GWTopts{1}.threshold0 = 0.5;
GWTopts{1}.precision  = 1e-4;
GWTopts{2}=GWTopts{1};
GWTopts{2}.ManifoldDimension = 0;
GWTopts{2}.precision  = 0.05;

%% Go parallel
if matlabpool('size')==0,
    matlabpool
end;


%% Load data from MFA models: it takes too long to generate these models so we've saved them
fprintf('\n Loading precomputed MFA models...');
load DataForPaper/script_test_MFA_DP_generative_min
fprintf('done.');


%% Create the GWT and SVD models (this is fast enough)
fprintf('\n Generating GWT and SVD models...');

if true,
for i = 1:length(DataSet),
    %% Generate the data set
    fprintf('\n\n---- Data set %d --------',i);
    X{i} = GenerateDataSets( DataSet(i).name, DataSet(i).opts );
        
    %Xcentered = bsxfun(@minus,X,mean(X,1));
    fprintf('\n Computing SVD...');
    Timings.SVD(i) = cputime;
    [~,S{i},V{i}]  = svd(X{i},0);
    Timings.SVD(i) = cputime-Timings.SVD(i);
    fprintf('done.');
   
    %% Compute GWT and the transform of the data
    fprintf('\n Computing GWT...');
    Timings.GWT(i) = cputime;
    GWT{i}     = GMRA(X{i},GWTopts{i});
    Data{i}    = FGWT_trainingData(GWT{i},X{i});
    Timings.GWT(i) = cputime-Timings.GWT(i);
    
    %% Build a family of models, one per GWT scale, for the density
    % Allocate memory    
    Entropy         = zeros(max(GWT{i}.Scales),1);
    GWTCost         = zeros(max(GWT{i}.Scales),1);
    
    fprintf('\n Constructing GWT and SVD models...');
    for j = 1:max(GWT{i}.Scales),       
        % Estimate density at scale j
        fprintf('\n Model at scale %d...',j);
        TimingsGWT_DensEst(i,j)     = cputime;
        DensEst_GWT{i,j}            = GWT_EstimateDensityAtFixedScale( GWT{i}, Data{i}, j );
        TimingsGWT_DensEst(i,j)     = cputime-TimingsGWT_DensEst(i,j);
        
        % Now compute the cost of using j scales
        subtree_idxs{i,j}     = [DensEst_GWT{i,j}.cp_idx,get_ancestors(GWT{i}.cp,DensEst_GWT{i,j}.cp_idx)];
        GWTCost(i,j)          = computeDictionaryCost( GWT{i}, subtree_idxs{i,j}(:) )*(1+size(X,1)/size(X,2));
        
        % Retain SVD dimensions so that SVD cost is the same as GWTCost(j)
        lSVDidxthres(i,j)     = min([ceil(GWTCost(i,j)/(sum(size(X)))),size(V{i},2)]);
        X_SVDproj             = X{i}*V{i}(:,1:lSVDidxthres(i,j));
        % Generate SVD model
        TimingsSVD_DensEst(i,j)  = cputime;        
        DensEst_SVD{i,j}      = kde(X_SVDproj','rot');
        TimingsSVD_DensEst(i,j)  = cputime-TimingsSVD_DensEst(i,j);
        %Entropy_SVD(i,j)      = entropy(DensEst_SVD{i,j});
        fprintf('done.');
    end;
    fprintf('\n---- done data set %d --------',i);
end;
Timings.GWT_DensEst = TimingsGWT_DensEst;
Timings.SVD_DensEst = TimingsSVD_DensEst;
end;

%% Validation phase for GWT, MFA and SVD models
fprintf('\n------------------------------------------');
fprintf('Validating GWT, SVD, MFA models...');
for i = 1:length(DataSet),
    DataSet(i).opts.NumberOfPoints = DataSet(i).opts.NumberOfPoints;    
    X_tmp            = GenerateDataSets( DataSet(i).name, DataSet(i).opts );                         % True validation data
    lValNumberOfPoints = size(X_tmp,1);
        
    for j = 1:size(DensEst_GWT,2),
        if isempty(DensEst_GWT{i,j}), continue;end;
        fprintf('\n Validation of model (Data=%i,Scale=%i):',i,j);
        for r = 1:pValN,
            fprintf('%d,',r);
            X_val{r}            = GenerateDataSets( DataSet(i).name, DataSet(i).opts );                      % True validation data
            lValNumberOfPoints = size(X_val{r},1);
            if j<=size(DensEst_MFA,2),
                X_MFA{r}        = draw_mfa(DensEst_MFA{i,j}, lValNumberOfPoints)';                                % Data generated from MFA model
            end;
            X_GWT{r}            = GenerateFromMultiplePlaneDens( DensEst_GWT{i,j}, lValNumberOfPoints )';       % Data generated from GWT model
            X_SVD{r}            = sample( DensEst_SVD{i,j}, lValNumberOfPoints );                               % Data generated from SVD model
            X_SVD{r}            = V{i}(:,1:lSVDidxthres(i,j))*X_SVD{r};X_SVD{r} = X_SVD{r}';
                        
            if pMaxDimForNNs<size(X_val{r},2),
                [Q,R] = qr(randn(size(X_val{r},2),pMaxDimForNNs),0);
            else
                Q = eye(size(X_val{r},2));
            end;
            if j<=size(DensEst_MFA,2),
                [HD_MFAval{i,j}(r),dists_P_MFAval{i,j}(r,:),dists_Q_MFAval{i,j}(r,:),idxs_P_MFAval{i,j}(r,:),idxs_Q_MFAval{i,j}(r,:),medianHaussDist_MFAval{i,j}(r)]    = HausdorffDistance(X_val{r}*Q,X_MFA{r}*Q);   % Hausdorff distance for MFA model
            end;
            [HD_GWTval{i,j}(r),dists_P_GWTval{i,j}(r,:),dists_Q_GWTval{i,j}(r,:),idxs_P_GWTval{i,j}(r,:),idxs_Q_GWTval{i,j}(r,:),medianHaussDist_GWTval{i,j}(r)]        = HausdorffDistance(X_val{r}*Q,X_GWT{r}*Q);   % Hausdorff distance for GWT model
            [HD_SVDval{i,j}(r),dists_P_SVDval{i,j}(r,:),dists_Q_SVDval{i,j}(r,:),idxs_P_SVDval{i,j}(r,:),idxs_Q_SVDval{i,j}(r,:),medianHaussDist_SVDval{i,j}(r)]        = HausdorffDistance(X_val{r}*Q,X_SVD{r}*Q);   % Hausdorff distance for SVD model
        end;
        tHD_valval = zeros(pValN);tHD_MFAMFA = zeros(pValN);tHD_GWTGWT = zeros(pValN);tHD_SVDSVD = zeros(pValN);
        for r1 = 1:pValN,
            fprintf('\n\t Internal model validation %d,',r1);
            for r2 = r1+1:pValN,
                if pMaxDimForNNs<size(X_val{r},2),
                    [Q,R] = qr(randn(size(X_val{r},2),pMaxDimForNNs),0);
                end;
                [tHD_valval(r1,r2),~,~,~,~,tHD_median_valval(r1,r2)] = HausdorffDistance(X_val{r1}*Q,X_val{r2}*Q);
                if j<=size(DensEst_MFA,2),
                    [tHD_MFAMFA(r1,r2),~,~,~,~,tHD_median_MFAMFA(r1,r2)] = HausdorffDistance(X_MFA{r1}*Q,X_MFA{r2}*Q);
                end;
                [tHD_GWTGWT(r1,r2),~,~,~,~,tHD_median_GWTGWT(r1,r2)] = HausdorffDistance(X_GWT{r1}*Q,X_GWT{r2}*Q);
                [tHD_SVDSVD(r1,r2),~,~,~,~,tHD_median_SVDSVD(r1,r2)] = HausdorffDistance(X_SVD{r1}*Q,X_SVD{r2}*Q);
            end;
        end;
        HD_valval(i,j) = mean(tHD_valval(find(tHD_valval~=0)));
        HD_median_valval(i,j) = mean(tHD_median_valval(find(tHD_median_valval~=0)));
        if j<=size(DensEst_MFA,2),
            HD_MFAMFA(i,j) = mean(tHD_MFAMFA(find(tHD_MFAMFA~=0)));
        end;
        HD_GWTGWT(i,j) = mean(tHD_GWTGWT(find(tHD_GWTGWT~=0)));
        HD_SVDSVD(i,j) = mean(tHD_SVDSVD(find(tHD_SVDSVD~=0)));
        HD_median_GWTGWT(i,j) = mean(tHD_median_GWTGWT(find(tHD_median_GWTGWT~=0)));
        HD_median_SVDSVD(i,j) = mean(tHD_median_SVDSVD(find(tHD_median_SVDSVD~=0)));
        fprintf('done.');
    end;
end;

fprintf('\n');

return;



%% Produce Figures for the paper
pFiguresPerDensityEstimatorX = 4;
pFiguresPerDensityEstimatorY = 8;
pFiguresPerDensityEstimator = pFiguresPerDensityEstimatorX*pFiguresPerDensityEstimatorY;
pI = 2;
pJ = 6;

lFigPerLine = pFiguresPerDensityEstimatorX;

lFig = zeros(28*lFigPerLine,28*lFigPerLine*3);
figure;
for l = 1:3,
    for k = 0:pFiguresPerDensityEstimator-1,        
        lFigIdxsX = 28*(fix(k/lFigPerLine))+1:28*(fix(k/lFigPerLine)+1);
        lFigIdxsY = 28*lFigPerLine*(l-1)+1+28*(rem(k,lFigPerLine))+1+l-1:28*lFigPerLine*(l-1)+1+28*(rem(k,lFigPerLine)+1)+l-1;        
        if l==1,
            z = GenerateFromMultiplePlaneDens( DensEst_GWT{pI,pJ}, 1,struct('alpha',0.75) );
            lFig(lFigIdxsX,lFigIdxsY) = reshape(z,[28,28]);
        elseif l==2,
            z = sample( DensEst_SVD{pI,pJ}, 1); z=V{i}(:,1:lSVDidxthres(pI,pJ))*z;z = z';
            lFig(lFigIdxsX,lFigIdxsY) = reshape(z,[28,28]);
        else
            z = draw_mfa(DensEst_MFA{pI,1},1);
            lFig(lFigIdxsX,lFigIdxsY) = reshape(z,[28,28]);
        end;
    end;   
end;
imagesc(lFig);hold on;    
colormap(map2); balanceColor;
for k=1:2,
    l=line([28*lFigPerLine*k,28*lFigPerLine*k],[1,size(lFig,2)]);set(l,'Color','k');
end;
set(gca,'xtick',[],'ytick',[]);
set(gca,'xtick',[28*lFigPerLine/2,28*lFigPerLine/2*3,28*lFigPerLine/2*5]);
set(gca,'xticklabel',{'GWT','SVD','MFA'});


%% Plot some more figures for the data sets: Hausdorff Distance
for i = 1:length(DataSet),
    for j = 1:size(DensEst_GWT,2),
        if ~isempty(DensEst_GWT{i,j}), lMaxScaleGWTModel = j; else break; end;
    end;
    figure;
    plot([mean(HD_valval(i,2))*ones(1,lMaxScaleGWTModel);cellfun(@mean,HD_GWTval(i,1:lMaxScaleGWTModel));cellfun(@mean,HD_SVDval(i,1:lMaxScaleGWTModel));cellfun(@mean,HD_MFAval(i,2))*ones(1,lMaxScaleGWTModel)]','--x');
    xlabel('Scale of the GWT model');ylabel('Hausdorff distances');
    hold on;
    plot([mean(HD_valval(i,1))*ones(1,lMaxScaleGWTModel);HD_GWTGWT(i,1:lMaxScaleGWTModel);HD_SVDSVD(i,1:lMaxScaleGWTModel);HD_MFAMFA(i,2)*ones(1,lMaxScaleGWTModel)]','-o');
    xlabel('Scale of the GWT model');ylabel('Hausdorff distances');
    if i==length(DataSet),legend({'Validation: dist.to training','GWT: dist. to training','SVD: dist. to training','MFA: dist. to training','Validation variability','GWT variability','SVD variability','MFA variability'},'Location','best');end;
    %title(DataSet(i).name);
end;



