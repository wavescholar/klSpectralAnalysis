function node = find_nearest_leaf_node(gW, x)

% This function finds the leaf node of the tree gW
% whose center is the closest to the point x (1 by D vector).
% Only two fields in the structure gW are actually needed: 
% (1) .cp (tree structure) and 
% (2) .Centers (cell array of the net centers)

[~, I] = min(sum((cat(1, gW.Centers{gW.LeafNodes}) - repmat(x, numel(gW.LeafNodes), 1)).^2, 2));

node = gW.LeafNodes(I);