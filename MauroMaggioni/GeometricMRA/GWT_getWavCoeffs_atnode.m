function [WavCoeffs,Wav_cp_idx,k_idxs,j_idxs] = GWT_getWavCoeffs_atnode( GWT, Data, cp_idx )

WavCoeffs   = [];
Wav_cp_idx  = [];
k_idxs      = [];
j_idxs      = [];

% Find the wavelets coefficients of points belonging to each child for the node cp_idx
Wav_cp_idx = find(GWT.cp==cp_idx);                                                   % First of all find the indices of the children

for k = 1:length(Wav_cp_idx),
    [k_idxs,j_idxs] = find(Data.Cel_cpidx==Wav_cp_idx(k));
    WavCoeffs{k} = [];
    for i = 1:length(k_idxs),
        WavCoeffs{k} = [WavCoeffs{k},Data.CelWavCoeffs{k_idxs(i),j_idxs(i)}'];
    end;
end;

return;