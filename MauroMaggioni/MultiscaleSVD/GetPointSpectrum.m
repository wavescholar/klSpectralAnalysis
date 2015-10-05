function [vStats,MSVD_Stats] = GetPointSpectrum( X, PtIdxs, NetsOpts )

%% Find the neighbors of the specified points
lAtria = nn_prepare(X);
[count,neighbors] = range_search(X,lAtria,PtIdxs,max(NetsOpts.Delta),-1);

%% Compute PCA for the specified points, at all the specified scales
MSVD_Stats.S = zeros(length(NetsOpts.Delta),length(PtIdxs),NetsOpts.SVDDim);
for k = 1:length(PtIdxs),
    for j = 1:length(NetsOpts.Delta),
        idxs = neighbors{k,1}(find(neighbors{k,2}<NetsOpts.Delta(j)));
        Xlocal = X(idxs,:);
        S = svd(bsxfun(@minus,Xlocal,mean(Xlocal,1)),0)/sqrt(length(idxs));
        MSVD_Stats.S(j,k,1:length(S)) = S;       
        MSVD_Stats.count(j,k) = length(idxs);
    end;
    lStats = EstimateDimFromSpectra( NetsOpts.Delta', MSVD_Stats.count(:,k), struct(), squeeze(MSVD_Stats.S(:,k,:)) );
    vStats.GoodScales{k} = lStats.GoodScales;
    vStats.DimEst(k)     = lStats.DimEst;
end;

MSVD_Stats.S_lbl = 'Sing. vals.';

return;