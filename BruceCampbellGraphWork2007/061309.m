    [x,y,z] = peaks(101);  
    plot3k({x y z},gradient(z),[-0.5 0.5],{'o',2},11,{'Peaks'},...
          'FontName','Arial','FontSize',18,'FontWeight','Bold');


% Move to motion on a two dimensional lattice, the code is much the same,
% all that we need to change is the dimension parameter
clear all 
N = 2;
M = 100000;
X = zeros(M,N);
for ii = 2:M
X(ii,:) = X(ii-1,:)+RandDir(N)';
end
hold on; 
%comet(X(:,1),X(:,2));
x =X(:,1);
y =X(:,2);
TRI = delaunay(x,y);
% Z = linkage(X,'average'); 
% %X = rand(1000,2);
% Y = pdist(X,'cityblock'); 
% Z = linkage(Y,'average'); 
% [H,T] = dendrogram(Z,'colorthreshold','default'); 
% set(H,'LineWidth',2) 
% T = clusterdata(X,'maxclust',10); 
% s = silhouette(X,T,'sqeuclid'); 
voronoi(x,y);
%plot(x,y,'.');
;
[vx, vy] = voronoi(x,y,TRI);
plot(vx(1,:),vy(2,:));
plot(x,y,'og'); 
k = convhull(x,y);
plot(x(k),y(k),'--rs','LineWidth',2,... 
'MarkerEdgeColor','k',... 
'MarkerFaceColor','g',... 
'MarkerSize',10) 
[idx, C, sumD, D] = kmeans(X, 10);
hold on; 
TRIkm = delaunay(C(:,1),C(:,2)); 
% 
triplot(TRIkm,C(:,1),C(:,2),'--s','LineWidth',2,... 
'MarkerEdgeColor','k',... 
'MarkerFaceColor','g',... 
'MarkerSize',12); 
plot(C(:,1),C(:,2),'r.','MarkerSize',30) 
voronoi(C(:,1),C(:,2),'r+'); 
triplot(TRI,x,y);
plot(X(:,1),X(:,2),'-g'); 
[vx, vy] = voronoi(C(:,1),C(:,2),'+b') 
voronoi(C(:,1),C(:,2),'+y'); 
plot(vx(1,:),vy(2,:));
% place random points in the plane 
C=pdist(X,'euclidean'); % compute distance between points 
C=exp(-.1*squareform(C)); 
% edge cost is a negative exponential of distance 
k=30; 
% # of partitions 
[ndx,Pi,cost]= grPartition(C,k,30);
colors=hsv(k); 
% plots points with appropriate colors 
colormap(colors)
cla
line(X(:,1),X(:,2),'MarkerEdgeColor',[0,0,0],'linestyle','none','marker','.'); 
for i=1:k 
line(X(find(ndx==i),1),X(find(ndx==i),2),... 'MarkerEdgeColor',colors(i,:),'linestyle','none','MarkerSize',22,'marker','.'); 
end
title(sprintf('Cost %g',cost)) 
colorbar 
hold  off 
K = convhull(x,y)
plot(x,y,'r+',vx,vy,'b-'),axis([0 1 0 1]) 
xx = -1:.05:1; yy = abs(sqrt(xx));
[x,y] = pol2cart(xx,yy);
'MarkerSize'
,10) 
,
'r-',x,y,'b+') 
%% Two dimensions, many drunkards
%
% It's kicking out time at the pub, where will our 50 drunkards end up? Or
% put more clinically let's simulate and plot many random walks leaving
% the same point.
K = 50; 
% Number of paths to simulate 
N = 2; 
% number of dimensions 
M = 10000; 
% number of steps 
X = zeros(M,N,K);
for
ii = 1:K 
for jj = 2:M 
X(jj,:,ii) = X(jj-1,:,ii) + RandDir(N)';
end 
end
figure;
axes; grid 
on; 
shg
for
ii = 1:K 
line(X(:,1,ii),X(:,2,ii),
'color',rand(1,3)); 
drawnow
pause(.05);
end
TRI = delaunay(x,y)
%% Drunkards... with Jet-packs
% AKA Random walks in 3 dimensions
clear 
all; 
hold 
on 
N = 3;
M = 100000;
X = zeros(M,N);
for
ii = 2:M 
X(ii,:) = X(ii-1,:)+RandDir(N)';
end
comet3(X(:,1),X(:,2),X(:,3));
plot3(X(:,1),X(:,2),X(:,3));
grid 
on; 
hold 
on 
data.X=X
%data= clust_normalize(data,'range');
plot3(data.X(:,1),data.X(:,2),data.X(:,2),'.') 
param.c=3;
param.vis=1;
param.val=2;
[idx, C, sumD, D] = kmeans(data.X, 3);
plot3(C(:,1),C(:,2),C(:,2),'go'); 
[idx, C, sumD, D] = kmeans(data.X, 4);
plot3(C(:,1),C(:,2),C(:,2),'ro'); 
[idx, C, sumD, D] = kmeans(data.X, 5);
plot3(C(:,1),C(:,2),C(:,2),'bo'); 
%[IDX, C] = KMEANS(X, K) returns the K cluster centroid locations in
% the K-by-P matrix C.
%
% [IDX, C, SUMD] = KMEANS(X, K) returns the within-cluster sums of
% point-to-centroid distances in the 1-by-K vector sumD.
%
% [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
% to every centroid in the N-by-K matrix D.
%
%result=kmeans(data,param);
result=kmeans(data,3);
hold 
on 
plot(result.cluster.v(:,1),result.cluster.v(:,2),
'ro') 
%%
% Multiple paths
K = 50; 
% Number of paths to simulate 
N = 3; 
% number of dimensions 
M = 1000; 
% number of steps 
X = zeros(M,N,K);
for
ii = 1:K 
for jj = 2:M 
X(jj,:,ii) = X(jj-1,:,ii) + RandDir(N)';
end 
end
figure;
axes; grid 
on; view(3); 
shg
for
ii = 1:K 
line(X(:,1,ii),X(:,2,ii),X(:,3,ii),'color',rand(1,3)); 
drawnow
pause(.05);
end