function idxs = get_partition_at_scale( GWT, j )

%
% Returns a cover of the set of leaves made by nodes at scale j, 
% plus any additional leaf that may be needed at scale coarser than j
% 

idxs = sort([find(GWT.Scales == j) GWT.LeafNodes(GWT.Scales(GWT.LeafNodes)<j)], 'ascend');

return;