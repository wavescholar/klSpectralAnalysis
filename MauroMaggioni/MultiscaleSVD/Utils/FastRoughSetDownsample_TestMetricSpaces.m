% This script is to test FastRoughSetDownsample in the case the data is an arbitrary metric space
% We start with data sets which are in fact Euclidean, to check that everything makes in this case.

clear lTime;

global DistInfo;

% First tests: create a circle with increasing number of points, size of nets independent of number of points
lN = [1000,5000];
fprintf('\n\n Increasing number of points, fixed number of points in net.');
for k = 1:length(lN),
    fprintf('\n Running data set of size %d...',lN(k));
    lParam=linspace(0,2*pi,lN(k));
    lScaleNearestPoints = 2*pi/lN(k);
    cX=[sin(lParam)',cos(lParam)']+0.2*lScaleNearestPoints*randn(lN(k),2);
    lOpts = struct('Delta',2*pi/30,'DoSVD',true,'AvoidSmallClusters',0,'NormalizedSVD',0,'DoD',1);
    % This runs the usual Euclidean version
    tic;vNet_1(k) = FastRoughSetDownsample( cX, lOpts );lTime(k)=toc;
    DisplayNetWithSVD( cX, vNet_1(k).Idxs, vNet_1(k).S, vNet_1(k).V );axis equal;drawnow;
    fprintf('%d pts in net, %.2f seconds.',length(vNet_1(k).Idxs),lTime(k));
    % Now run the metric space version
    [DistInfo.count,DistInfo.idxs, DistInfo.dists] = nrsearch( cX', cX', 100 );
    AddSortedIdxsToDistInfo
    lOpts.DistInfo = DistInfo; lOpts.SVDDim = 10;
    tic;vNet_metric(k) = FastRoughSetDownsample( cX, lOpts );lTime_metric(k)=toc;
    fprintf('\n Metric space version: %d pts in net, %.2f seconds.',length(vNet_metric(k).Idxs),lTime_metric(k));
end;

fprintf('\n');
