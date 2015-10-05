function W = get_wavCoeffs_node(gW, Data, n)

% get the wavelet coefficients of all points in the node n

[~, offspringLeafNodes] = get_offspring(gW.cp, n);

imap(gW.LeafNodes) = 1:length(gW.LeafNodes);
W = cat(1, Data.CelWavCoeffs{imap(offspringLeafNodes), gW.Scales(n)});
