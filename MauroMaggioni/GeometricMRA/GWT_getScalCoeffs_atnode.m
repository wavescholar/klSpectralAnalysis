function [ScalCoeffs,k_idxs,j_idxs] = GWT_getScalCoeffs_atnode( GWT, Data, cp_idx )

% Find the scaling coefficients of points belonging to each node at scale j
[k_idxs,j_idxs] = find(Data.Cel_cpidx==cp_idx);
ScalCoeffs = [];

for i = 1:length(k_idxs),
    ScalCoeffs = [ScalCoeffs,Data.CelScalCoeffs{k_idxs(i),j_idxs(i)}'];
end;

return;