%%
% Noisy circle example
%
lN = 1000;
lParam=linspace(0,2*pi,lN);cX=[sin(lParam)',cos(lParam)']+0.01*randn(lN,2);cX=cX';
lOpts = struct('Delta',0.1,'DoSVD',true,'AvoidSmallClusters',1,'FastNNSearcher','covertree');
profile on;tic;
vNet = FastRoughSetDownsample( cX, lOpts );
toc;profile viewer;
DisplayNetWithSVD( cX, vNet.Idxs, vNet.S, vNet.V );axis equal;gca; hold on;plot(cX(1,find(vNet.InvIdxs==10)),cX(2,find(vNet.InvIdxs==10)),'k.');plot(cX(1,vNet.Idxs(10)),cX(2,vNet.Idxs(10)),'ro');
figure;plot(cX(1,:),cX(2,:),'c.');hold on;plot(cX(1,find(vNet.InvIdxs==10)),cX(2,find(vNet.InvIdxs==10)),'k.');hold on;plot(cX(1,vNet.Idxs(10)),cX(2,vNet.Idxs(10)),'ro');



%%
%
% SPEED TEST
%
%

clear lTime;

% First tests: create a circle with increasing number of points, size of nets independent of number of points
lN = [1000,5000,10000,50000,100000];
fprintf('\n\n Increasing number of points, fixed number of points in net.');
for k = 1:length(lN),
    fprintf('\n Running data set of size %d...',lN(k));
    lParam=linspace(0,2*pi,lN(k));
    lScaleNearestPoints = 2*pi/lN(k);
    cX=[sin(lParam)',cos(lParam)']+0.2*lScaleNearestPoints*randn(lN(k),2);
    lOpts = struct('Delta',2*pi/300,'DoSVD',true,'AvoidSmallClusters',0,'NormalizedSVD',0,'DoD',1);
    tic;vNet_1(k) = FastRoughSetDownsample( cX, lOpts );lTime(k)=toc;
    DisplayNetWithSVD( cX, vNet_1(k).Idxs, vNet_1(k).S, vNet_1(k).V );axis equal;drawnow;
    fprintf('%d pts in net, %.2f seconds.',length(vNet_1(k).Idxs),lTime(k));
end;

figure;plot(lN,lTime);title('Computation time as function of size of data set');

fprintf('\n');

clear lTime;

% Create a circle with increasing number of points, size of nets dependent of number of points
fprintf('\n\n Increasing number of points, number of points in net also increasing (proportionally).');
for k = 1:length(lN),
    fprintf('\n Running data set of size %d...',lN(k));
    lParam=linspace(0,2*pi,lN(k));
    lScaleNearestPoints = 2*pi/lN(k);
    cX=[sin(lParam)',cos(lParam)']+0.2*lScaleNearestPoints*randn(lN(k),2);
    lOpts = struct('Delta',lScaleNearestPoints *100,'DoSVD',true,'AvoidSmallClusters',0,'NormalizedSVD',0);
    tic;vNet_2 = FastRoughSetDownsample( cX, lOpts );lTime(k)=toc;
    DisplayNetWithSVD( cX, vNet_2.Idxs, vNet_2.S, vNet_2.V );axis equal;drawnow;
    fprintf('%d pts in net, %.2f seconds.',length(vNet_2.Idxs),lTime(k));
end;

figure;plot(lN,lTime);title('Computation time as function of size of data set');

fprintf('\n');

clear lTime;
% Now we test increases in extrinsic dimensionality
lD = [2,4,6,8,10,15,20,40,60,100,200,500,1000];
fprintf('\n\n Increasing dimensions, fixed number of points.');
for k = 1:length(lD),
    fprintf('\n Running data set in dimension %d...',lD(k));
    % Compute a random rotation in lD dimensions
    lT = randn(lD(k),lD(k)); [lT,lR] = qr(lT);
    lParam=linspace(0,2*pi,lN(1));
    lScaleNearestPoints = 2*pi/lN(1);
    cX=[sin(lParam)',cos(lParam)']+0.2*lScaleNearestPoints*randn(lN(1),2);
    cX=[cX,zeros(size(cX,1),lD(k)-2)]*lT;
    lOpts = struct('Delta',2*pi/300,'DoSVD',true,'AvoidSmallClusters',0,'NormalizedSVD',0,'SVDDim',3);
    tic;vNet_3 = FastRoughSetDownsample( cX, lOpts );lTime(k)=toc;
    DisplayNetWithSVD( cX, vNet_3.Idxs, vNet_3.S, vNet_3.V );axis equal;drawnow;
    fprintf('%d pts in net, %.2f seconds.',length(vNet_3.Idxs),lTime(k));
end;

figure;plot(lD,lTime);title('Computation time as function of extrinsic dimension');

fprintf('\n');

clear lTime;
% Now we test increases in extrinsic dimensionality
lD = [2,4,6,8,10,15,20,40,60,100,200,500,1000];
fprintf('\n\n Increasing dimensions, increase number of points proportionally to dimension.');
for k = 1:length(lD),
    fprintf('\n Running data set in dimension %d...',lD(k));
    % Compute a random rotation in lD dimensions
    lT = randn(lD(k),lD(k)); [lT,lR] = qr(lT);
    lNpts = round(lN(1)*lD(k)/10);
    lParam=linspace(0,2*pi,lNpts);
    lScaleNearestPoints = 2*pi/lNpts;
    cX=[sin(lParam)',cos(lParam)']+0.2*lScaleNearestPoints*randn(lNpts,2);
    cX=[cX,zeros(size(cX,1),lD(k)-2)]*lT;
    lOpts = struct('Delta',2*pi/300,'DoSVD',true,'AvoidSmallClusters',0,'NormalizedSVD',0,'SVDDim',3);
    tic;vNet_3 = FastRoughSetDownsample( cX, lOpts );lTime(k)=toc;
    DisplayNetWithSVD( cX, vNet_3.Idxs, vNet_3.S, vNet_3.V );axis equal;drawnow;
    fprintf('%d pts in net, %.2f seconds.',length(vNet_3.Idxs),lTime(k));
end;

figure;plot(lD,lTime);title('Computation time as function of extrinsic dimension d, and linear number of points d');

fprintf('\n');

return