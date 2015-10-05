%% This scripts tests simple modeling of densities, comparing models in svd space vs. GWT space

clear all;%close all;clc; clear global

pN      = 4000;
pValN   = pN;
pVarN   = 20;
pTrials = 1;

%% Generate Data set, and a validation data set
XName = 'BMark_MNIST'; XNickName = 'MNIST Digits';
XOpts = struct('NumberOfPoints',pN,'MnistOpts',struct('Sampling', 'RandN', 'QueryDigits',7, 'ReturnForm', 'vector'));


%% Go parallel
if matlabpool('size')==0,
    matlabpool
end;

%% Generate the data set                
X = GenerateDataSets( XName, XOpts);

for p = 1:pTrials,
    Xval = GenerateDataSets( XName, XOpts);
    HausDistXXval(p) = HausdorffDistance(X,Xval);
end;
HausDistXXval = mean(HausDistXXval);

Xcentered = bsxfun(@minus,X,mean(X,1));
Timings.SVD = cputime;
[U,S,V] = svd(Xcentered,0);
Timings.SVD = cputime-Timings.SVD;

%% Compute GWT
Timings.GWT = cputime;
GWTopts = struct('GWTversion',0);
GWTopts.ManifoldDimension = 0;
GWTopts.threshold1 = 1e-3; 
GWTopts.threshold2 = .1; 
GWTopts.addTangentialCorrections = false;
GWTopts.sparsifying = false;
GWTopts.splitting = false;
GWTopts.knn = 30;
GWTopts.knnAutotune = 20;
GWTopts.smallestMetisNet = 10;
GWTopts.verbose = 1;
GWTopts.shrinkage = 'hard';
GWTopts.avoidLeafnodePhi = false;
GWTopts.mergePsiCapIntoPhi  = true;
GWTopts.coeffs_threshold = 0;
GWTopts.errorType = 'absolute';
GWTopts.threshold0 = 0.5; 
GWTopts.precision  = 1e-2;

GWT     = GMRA(X,GWTopts);
Data    = FGWT_trainingData(GWT,X);
Timings.GWT = cputime-Timings.GWT;

%% Build a family of models, one per GWT scale, for the density
% Allocate memory
HausDist        = zeros(max(GWT.Scales),pTrials);
Entropy         = zeros(max(GWT.Scales),1);
GWTCost         = zeros(max(GWT.Scales),1);

for j = 1:max(GWT.Scales),
    fprintf('\n Generating model at scale %d...',j);
    % Estimate density at scale j
    Timings.GWT_EstDens(j) = cputime;
    [DensEst{j},cp_idx{j}] = GWT_EstimateDensityAtFixedScale( GWT, Data, j );
    Timings.GWT_EstDens(j) = cputime-Timings.GWT_EstDens(j);
    for p = 1:pTrials,
        % Generate new data
        Xnew{j} = GenerateFromMultiplePlaneDens( DensEst{j}, pValN )';        % MM_DBG: figure;plot3(Xnew{j}(:,1),Xnew{j}(:,2),Xnew{j}(:,3),'.');hold on;p=plot3(GWT.X(:,1),GWT.X(:,2),GWT.X(:,3),'r.');set(p,'MarkerSize',0.5);
        % Evaluate a distance between the new samples and the validation samples
        HausDist(j,p) = HausdorffDistance(Xnew{j},Xval);        
    end;
    % Estimate entropy
    Entropy(j)          = DensEst{j}.entropy;
    % Now compute the cost of using j scales
    subtree_idxs{j}     = [cp_idx{j},get_ancestors(GWT.cp,cp_idx{j})];
    GWTCost(j)          = computeDictionaryCost( GWT, subtree_idxs{j}(:) )*(1+size(X,1)/size(X,2));
    
    % Retain SVD dimensions so that SVD cost is the same as GWTCost(j)
    lSVDidxthres(j) = min([ceil(GWTCost(j)/(sum(size(X)))),size(V,2)]);
    X_SVDproj       = Xcentered*V(:,1:lSVDidxthres(j));
    Timings.SVD_Dens(j) = cputime;
    DensEstSVD{j}       = kde(X_SVDproj','rot');
    Timings.SVD_Dens(j) = cputime-Timings.SVD_Dens(j);
    EntropySVD(j)       = entropy(DensEstSVD{j});
    XSVDnew{j}          = sample(DensEstSVD{j},pValN);
    XSVDnew{j}          = V(:,1:lSVDidxthres(j))*XSVDnew{j};
    XSVDnew{j}          = XSVDnew{j}';
    HausDistSVD(j)      = HausdorffDistance(XSVDnew{j},Xval);
    
    fprintf('done.');
end;


%% Now estimate variability of the model in terms of Hausdorff distances of genereated data sets
% Generate several new data sets
HausDistGeneratedData   = zeros(length(DensEst),pVarN,pVarN);
HausDistToVal           = zeros(length(DensEst),pVarN);
for j = 1:length(DensEst),
    for i = 1:pVarN;
        Z{i}=GenerateFromMultiplePlaneDens( DensEst{j}, pValN )';
    end;
    for i=1:pVarN;
        for l=i+1:pVarN;
            HausDistGeneratedData(j,i,l)=HausdorffDistance(Z{i},Z{l});
        end;
        HausDistToVal(j,i)=HausdorffDistance(Z{i},Xval);
    end;
end;

HausDistGeneratedDataSVD   = zeros(length(DensEst),pVarN,pVarN);
HausDistToValSVD           = zeros(length(DensEst),pVarN);
for j = 1:length(DensEst),
    for i = 1:pVarN;
        Z{i}=sample( DensEstSVD{j}, pValN )';
        Z{j}      = V(:,1:lSVDidxthres(j))*Z{j};
        Z{j}      = Z{j}';
    end;
    for i=1:pVarN;
        for l=i+1:pVarN;
            HausDistGeneratedDataSVD(j,i,l)=HausdorffDistance(Z{i},Z{l});
        end;
        HausDistToValSVD(j,i)=HausdorffDistance(Z{i},Xval);
    end;
end;


%% Plot some figures
figure;plot(mean(HausDist,2)/HausDistXXval);hold on;plot(HausDistSVD/HausDistXXval,'r');legend({'GWT','SVD'});title('Normalized Hausdorff distance between generated points and validation points');xlabel('j');
figure;plot(Entropy);hold on;plot(EntropySVD,'r');legend({'GWT','SVD'});title('Entropy of estimated distributions');xlabel('j');

figure;imagesc(reshape(Xnew{6}(2,:),[28,28]));colormap(gray);title('Example of reconstructed digit from GWT model');
figure;imagesc(reshape(XSVDnew{6}(2,:),[28,28]));colormap(gray);title('Example of reconstructed digit from SVD model');

%%
return;