function p = EvalFromMultiplePlaneDens( GWT,DensEst,Y )

% Computes density at all points in Y (each row is a data point) at all scales of the GWT and DensEst
%
% Example:
% for j = 1:max(GWT.Scales),[DensEst{j},cp_idx{j}] = GWT_EstimateDensityAtFixedScale( GWT, Data, j );end;
% p = EvalFromMultiplePlaneDens( GWT,DensEst,GWT.X );figure;imagesc(p);colorbar;sum(log(p),1)
%

p = zeros(size(Y,1),max(GWT.Scales));

for k = 1:size(Y,1),
    yhat = FGWT(GWT,Y(k,:));
    for j = 1:length(DensEst),
        if isempty(DensEst{j}), continue; end;
        [idx,idx_chain,idx_dens] = intersect(yhat.chain,DensEst{j}.cp_idx);
        p(k,j) = evaluate(DensEst{j}.Density(idx_dens),yhat.ScalCoeffs{idx_chain}');    %DensEst{j}.pi(idx_dens)*evaluate(DensEst{j}.Density(idx_dens),yhat.ScalCoeffs{idx_chain}');
    end;
end;

return;