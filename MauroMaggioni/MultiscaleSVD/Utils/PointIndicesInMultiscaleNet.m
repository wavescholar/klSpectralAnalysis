function vPointIdxs = PointIndicesInMultiscaleNet( cPointIdxs, cNets )

% OUT:
%   vPointIdxs(i,j) is the cluster index in the net at scale j of cPointIdxs(i)

vPointIdxs(length(cPointIdxs),length(cNets)) = uint32(0);

for lj = 1:length(cNets),
    vPointIdxs(:,lj) = uint32(cNets(lj).AbsInvIdxs(cPointIdxs));
end;

return;