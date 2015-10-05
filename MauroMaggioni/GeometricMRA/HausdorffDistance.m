function [HaussDist,dists_P,dists_Q,idxs_P,idxs_Q,medianHaussDist] = HausdorffDistance( P,Q )

% P,Q : two point clouds (#pts * dim matrices)
% Expects matlabpool already open for parallelism


[idxs_P,dists_P] = mindistancessq( P, Q );  dists_P=sqrt(dists_P);
[idxs_Q,dists_Q] = mindistancessq( Q, P );  dists_Q=sqrt(dists_Q);

HaussDist       = max([max(dists_P),max(dists_Q)]);
medianHaussDist = mean([median(dists_P),median(dists_Q)]);

return;

function [idxs_P,dists_P] = mindistancessq( P, Q )

dists_P     = zeros(size(P,1),1);
idxs_P      = zeros(size(P,1),1);

parfor k = 1:size(P,1),
    Xtmp = Q;
    diffs = bsxfun(@minus,Q,P(k,:));
    norms = sum(diffs.^2,2);
    [dists_P(k),idxs_P(k)] = min(norms);
end;

return;

% % % % % If P (and therefore Q) have dimension 1, we need to fix a bug in nn_search
% % % % if size(P,2)==1,
% % % %     P = [P,zeros(size(P,1),1)];
% % % %     Q = [Q,zeros(size(P,1),1)];
% % % % end;
% % % % 
% % % % Atria_P = nn_prepare( P );
% % % % Atria_Q = nn_prepare( Q );
% % % % 
% % % % % Find nearest neighbor in P for every point in Q
% % % % [idxs_P,dists_P] = nn_search( P, Atria_P, Q, 1);
% % % % % Find nearest neighbor in P for every point in Q
% % % % [idxs_Q,dists_Q] = nn_search( Q, Atria_Q, P, 1);
% % % % 
% % % % HaussDist = max([max(dists_P),max(dists_Q)]);
% % % % 
% % % % return;
% % % % 
% % % % 
% % % % 
% % % % 
% % % % % % % [count,idxs,dists] = nrsearch(P, Q, 1, 0,[],struct('XIsTransposed',true,'ReturnAsArrays',true));
% % % % % % % [count,idxs,dists2] = nrsearch(Q, P, 1, 0,[],struct('XIsTransposed',true,'ReturnAsArrays',true));
% % % % % % % % [idxs, dists] = ANNsearch(cX, cXq, NN, opts.Tolerance);   
% % % % % % % 
% % % % % % % HaussDist2 = max(max(dists),max(dists2)),
% % % % 
% % % % 

