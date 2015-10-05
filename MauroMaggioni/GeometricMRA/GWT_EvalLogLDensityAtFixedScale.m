function [LogL,LogL_w,EvaluateDens,EvaluateDens_w] = GWT_EvalLogLDensityAtFixedScale( DensEst, GWT, Data )

%% Compute the LogL for the density of scaling function coefficients
for k = 1:length(DensEst.cp_idx),
    % Get the scaling function and wavelet coefficients
    lScalCoeffs = GWT_getScalCoeffs_atnode ( GWT, Data, DensEst.cp_idx(k) );
    lWavCoeffs  = GWT_getWavCoeffs_atnode  ( GWT, Data, DensEst.cp_idx(k) );
    % Compute the corresponding LogL
    LogL(k)     = evalAvgLogL(DensEst.Density(k),lScalCoeffs);
    EvaluateDens{k}= evaluate(DensEst.Density(k),lScalCoeffs); 
    for r = 1:length(lWavCoeffs),
        LogL_w(k,r)  = evalAvgLogL(DensEst.Density_w(k,r),lWavCoeffs{r});
        EvaluateDens_w{k,r} = evaluate(DensEst.Density_w(k,r),lWavCoeffs{r});
    end;
end;

return; 
