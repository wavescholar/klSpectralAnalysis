%%%%%%%%%%%%%%%%%%%%%%%%%%%CLUSTER In 2d  --
MU1 = [1.5 2.5];
SIGMA1 = [2 0; 0 .5];
MU2 = [-2 -3];
SIGMA2 = [1 0; 0 1];
MU3 = [0 0];
SIGMA3 = [1 0; 0 1];
X = [mvnrnd(MU1,SIGMA1,100);mvnrnd(MU2,SIGMA2,100);]
    mvnrnd(MU3,SIGMA3,100)];
scatter(X(:,1),X(:,2),30,'.')
hold on
options = statset('Display','final');
obj = gmdistribution.fit(X,3,'Options',options);
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
idx = cluster(obj,X);
cluster1 = X(idx == 1,:);
cluster2 = X(idx == 2,:);
cluster3 = X(idx == 3,:);
h1 = scatter(cluster1(:,1),cluster1(:,2),30,'r.');
h2 = scatter(cluster2(:,1),cluster2(:,2),30,'g.');
h3 = scatter(cluster3(:,1),cluster3(:,2),30,'b.');
legend([h1 h2 h3 ],'Cluster 1','Cluster 2','Cluster 3','Location','NW')
P = posterior(obj,X);
figure
scatter(X(:,1),X(:,2),30,P(:,1),'.')
hb = colorbar;
ylabel(hb,'Component 1 Probability')
D = mahal(obj,X);
figure
scatter(X(:,1),X(:,2),30,D(:,1),'.')
hb = colorbar;
ylabel(hb,'Mahalanobis Distance to Component 1');

[row_clust_idx, col_clust_idx,y_index,x_index]= SpectralCoClustering(X,2,1,['c1' 'c2']);



%%%%%%%%%%%%%%%%%%%%%END CLuster TB
%D=compute_distance_matrix(X);

D=pdist(X,'euclidean');
%PDIST Pairwise distance between observations.
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the N-by-P data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is a
%   1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
%   observations in X.
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       function      - A distance function specified using @, for
%                       example @DISTFUN

%drtoolbox
%mappedX = compute_mapping(X, 'Laplacian', no_dims, 7);	


%Make a gram matrix
G_11=gram(cluster1,cluster1,'gauss',1);
G_12=gram(cluster1,cluster2,'gauss',1);
G_13=gram(cluster1,cluster3,'gauss',1);

G_22=gram(cluster2,cluster2,'gauss',1);
G_23=gram(cluster2,cluster3,'gauss',1);
G_33=gram(cluster3,cluster3,'gauss',1);

%Linear kernel: G = gram(X1, X2, 'linear')
%           which is parameterless
% Gaussian kernel: G = gram(X1, X2, 'gauss', s)
%           where s is the variance of the used Gaussian function (default = 1).
% Polynomial kernel: G = gram(X1, X2, 'poly', R, d)
%           where R is the addition value and d the power number (default = 0 and 3)


G1=[G_11 G_12 G_13];
G2=[ G_12' G_22 G_23];
G3=[G_13' G_23' G_33];

G=[G1' G2' G3']';


% 
% GEV_l=G>0.2;% & G<.8;
 GEV_u=G>0.8 ;%& G<.8;
 GEV_u=GEV_u+GEV_u'
% %[{'square'} | 'heart' | 'sphere' | 'circle' | 'ballN'* | 'cubeN'*]
XY = gursoy_atun_layout(sparse(GEV_u*1.0),'topology','circle');
XY = kamada_kawai_spring_layout(sparse(GEV_u*1.0));
% XY=fruchterman_reingold_force_directed_layout(sparse(GEV_u*1.0));
% wgplot(GEV_l,X);
% AB=sparse(double(Net));
% XB = gursoy_atun_layout(AB,'topology','heart');gplot(AB,X,'.-');
wgplot(GEV_u,XY);










%creates the tree using the specified method. 
%Methods differ from one another in how they measure the distance between clusters. 
% %Available methods 
% 'average'	Unweighted average distance (UPGMA).
% 'centroid'	Centroid distance (UPGMC). Y must contain Euclidean distances.
% 'complete'	Furthest distance.
% 'median'	Weighted center of mass distance (WPGMC). Y must contain Euclidean distances.
% 'single'	Shortest distance. This is the default.
% 'ward'	Inner squared distance (minimum variance algorithm). Y must contain Euclidean distances.
% 'weighted' Weighted average distance (WPGMA).
Z = linkage(D,'average'); 

%X = rand(1000,2);
Y = pdist(X,'cityblock');
Z = linkage(Y,'average');
[H,T] = dendrogram(Z,'colorthreshold','default');
set(H,'LineWidth',2)

T = clusterdata(X,'maxclust',3); 

s = silhouette(X,T,'sqeuclid');

N=2048%4096;
%Generate a Winger - GOE
G = randn(N)/sqrt(N); 
A = (G+G')/sqrt(2); 
C1=eig(A);
C_Epdf1=hist(C1,200);


%Generate Wishart Matrix
G = randn(N,2*N)/sqrt(2*N); 
B = G*G';
C=eig(B);
EmpiricalPDF= hist(C,200);
plot(EmpiricalPDF);


%PCA Analysis
stdr = std(A);
sr = A./repmat(stdr,N,1);
[coefs,scores,variances,t2] = princomp(sr);
pca(A);
C=eig(A);


%http://www.stanford.edu/~dgleich/demos/matlab/spectral/spectral.html
clear all;
N=96;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = randperm(N);
gs = N/2;
G1 = x(1:gs);
G2 = x(gs+1:end);
% Decide on the probabilities of edges within each group and between the two groups. 
% Because we are planting a partition, the probabilities of edges between the groups should be much 
% lower than the probability of edges within each group. Suppose that group 1 is a little more tightly connected than group 2.
p_G1 = 0.23;
p_G2 = 0.43;
p_Inbetween = 0.04;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(N-gs,N-gs) < p_G2;
A(G1, G2) = rand(gs, N-gs) < p_Inbetween;
At = triu(A,1);

A = At -At';

L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);


[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(A);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);


[V D] = eigs(L2, 2, 'SA');
D(2,2);%ans = 46.7158
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:,2));






Dist = zeros(N,N);
for i=1:(N-1);
  for j=i+1:N;
    Dist(i,j) = norm( A(i,:) - A(j,:) );
  end;
end;
Dist = Dist + Dist';
threshold=9;
Ad = Dist ;%< threshold;
Ad = Ad - eye(N);
m = sum(sum(Ad))/2;

% find the incidence matrix
Ai = zeros(N,N);
l = 0;
for i=1:(N-1);
  for j=i+1:N;
    if Ad(i,j)>0.5
      l = l + 1;
      Ai(i,l) =  1;
      Ai(j,l) = -1;
    end;
  end;
end;
Ais = sparse(Ai);

spy(Ai);
%L = laplacian(A);    %BBCREVISIT - where did he get the laplacian
%from?
%L= DELSQ (A);
L1 = del2(A);%L=compute_laplacian(A);%
L2=laplacian(A);%
L3=compute_combinatorial_laplacian(A);
[V D] = eigs(L, 2, 'SA');
D(2,2);%ans = 46.7158
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:,2));
spy(A(p,p));draw_dot(A(p,p));
%[ignore p] = sort(V(:,2));
%spy(A(p,p));
%Let's do an MDS on the graph adjacency
 [points]=mds(A,2);%metric,iterations,learnrate)
 %D=all_shortest_paths(A);
 plot(points(:,1),points(:,2));
[ w_mh   ] = mh(A);
plotgraphBruce(A,points,w_mh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Generation of cluster grapg [my own term ]
% Three cluseters with a cut between each
clear all;
N=3*296;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.13;
p_G2 = 0.43;
p_G3 = 0.23;
p_Inbetween = 0.06;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;


[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));

B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween;
 A(G2, G3) = B_2 < p_Inbetween/3;
 A(G1, G3) = B_3 < p_Inbetween/5;

A = triu(A,1);
A = A + A';spy(A);
L = del2(A);spy(L);
[V D] = eigs(L, 2, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2)
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

[x y]=draw_dot(A(p,p));
gplot(A(p,p), [x' y'], '.-')

C_coeff=clustering_coefficients((sparse(A)));full(A)
hist(C_coeff,100);C_mu=mean(C_coeff);

Net = SFNG(N+20, 5, A);
C_coeff=clustering_coefficients(sparse(double(Net)));
hist(C_coeff,100);
C_mu=mean(C_coeff);

[x y]=draw_dot(Net);
gplot(Net, [x' y'], '.-');

complete_subgraphs=maximalCliques(B_1complete_subgraphs);

%%%%%%%%%%%%%%%%%FanChung Find
%http://www.math.ucsd.edu/~fan/algo/sigma.html
[W, L, sigma] = sigma_calc_weights(B_1, 0);
 sigma_write_dot(example_graph, W, 'example.dot', 0, 0); 



[ w_fdla ] = fdla(sparse(A(p,p)));
[ w_fmmc ] = fmmc(Net);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(Net);
[ w_mh   ] = mh(Net);

%PL_Equation = PLplot((Net));
%CNet(Net);
%Gerneating Erdos Reyni graphs
 clear all;
 N=50;
 A = rand(N+20,N+20);
 A = triu(A,1);
 p=log(N+20)/N+20-1/(N+20^2);%.02;%
 thres_conn=log(N+20)/N+20;
 A = A + A';
 G =  A < p;spy(A);draw_dot(A);
 thres_giant = 1/(N+20-1);
[x y] = draw_dot(G<thres_giant);
 gplot(A<thres_giant, [x' y'], '.-') 
[x y] = draw_dot(G<thres_conn);
 gplot(G<thres_conn, [x' y'], '.-')
C_coeff=clustering_coefficients(sparse(double(A)));spy(sparse(A));
hist(C_coeff,100);C_mu=mean(C_coeff);

n = 1500;
p = 0.04; 
G = rand(n,n) < p;
G = triu(G,1);
G = G + G';
thres_giant = 1/(n-1);
[x y] = draw_dot(G);
gplot(G, [x' y'], '.-'); 


%Now Const a random G from V -

% Given a square matrix M, the goal is to find a vector (with dii > 0)
% such that ||DMD^{-1}||_F is minimized, where D = diag(d).
% The problem can be cast as an unconstrained geometric program:
%           minimize sqrt( sum_{i,j=1}^{n} Mij^2*di^2/dj^2 )
%
% formulating the problem as a GP
cvx_setup
cvx_begin gp
  variable d(N)
  minimize( sqrt( sum( sum( diag(d.^2)*(A.^2)*diag(d.^-2) ) ) ) )
  % Alternate formulation: norm( diag(d)*abs(M)*diag(1./d), 'fro' )
cvx_end

A_min_spectral_radius_F=diag(d)*A* inv(diag(d));
L2=del2(A_min_spectral_radius_F);
;%d* ( 1 ./ d)';
[V2 D2] = eigs(L2, 2, 'SR');
D2(2,2);
plot(V2(:,2), '.-');
plot(sort(V2(:,2)), '.-');spy(A_min_spectral_radius_F);
[ignore p] = sort(V2(:,2));
spy(A_min_spectral_radius_F(p,p));draw_dot(A_min_spectral_radius_F(p,p));

hold on;
% displaying results
D = diag(d);
disp('The matrix D that minimizes ||DMD^{-1}||_F is: ');
disp(D);
disp('The minimium Frobenius norm achieved is: ');
disp(norm(D*M*inv(D),'fro'));
disp('while the Frobunius norm of the original matrix M is: ');
disp(norm(M,'fro'));







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BGL %%%%%%%%%%%%%%%%CORE NUM

% For fun, let's see how many cores there are in a road network.  Vertices
% in a 1-core in a road network have at least one path between them
% (assuming the underlying network is connected).  Vertices in a 2-core
% have at least two paths between them.  
load('graphs/minnesota.mat');A=spones(A); % convert to unweighted
Net = SFNG(2642+100, 2, A);
[cn csz]=core_numbers(A);  cs=unique(cn);spy(A);
 %[A,xy] = VxE_to_AdjacencyMatrix_and_XY_points(V,E);
 X = gursoy_atun_layout(A,'topology','heart');
 %X = kamada_kawai_spring_layout(A);%fruchterman_reingold_force_directed_layout(A);
gplot(A,X,'.-');
AB=sparse(double(Net));
XB = gursoy_atun_layout(AB,'topology','heart');gplot(AB,X,'.-');

C_coeff=clustering_coefficients(AB)
hist(C_coeff,100);

%%%%%%%%%%%%%%%%%%%%%%%%%Plot Cores
clf; hold on; wh=ones(1,3); colors={0.85*wh,0.7*wh,0.55*wh,0.4*wh};
for c=unique(cn')
  m=cn>=c; cl=colors{c+1}; lw=16*1.5^(max(cn)-c);
  xym=xy(m,:);is=convhull(xym(:,1),xym(:,2));
  h=fill(xym(is,1),xym(is,2),cl);set(h,'LineWidth',lw,'EdgeColor',cl);
end
gplot(A,X,'k-');  plot(X(:,1), X(:,2),'r.','MarkerSize',24);
text(X(:,1)+0.1, X(:,2)+0.1, num2str((1:21)'));set(gcf,'Color',[1,1,1]);
xlim([-1,10]);ylim([-2,7]);axis off; hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End plot cors




gplot(G,X,'.-');%draw_dot(full(A));
arrayfun(@(v,c) fprintf('core(%2i) = %4i\n',v,c),cs,csz(cs+1));
% Now, you might complain that there are certain vertices in this graph
% that simply chain a path in the road network so it draws correctly.  The
% following code removes all vertices of degree two and connects the
% end-points of the degree two vertex directly.  It applies this prodecure
% iteratively until all the vertices of degree two are gone.  At the end,
% we compute the core nubmers again.
d2v=find(sum(A,2)==2);
while ~isempty(d2v)
    for v=d2v'
        l=find(A(:,v));A(l(1),v)=0;A(l(2),v)=0;A(v,l(1))=0;A(v,l(2))=0;
        A(l(1),l(2))=1;A(l(2),l(1))=1;
    end
    d2v=find(sum(A,2)==2);
end
max(core_numbers(A))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END BGL 

%http://www.stanford.edu/~dgleich/demos/matlab/random_graphs/erdosreyni.html
%Gerneating Erdos Reyni graphs
clear all
n = 200;
A = rand(n,n);
A = triu(A,1);
p=log(n)/n-1/(n^2);%.02;%
thres_conn=log(n)/n;
A = A + A';
G =  A < p;
thres_giant = 1/(n-1);
[x y] = draw_dot(G<thres_giant);
gplot(A<thres_giant, [x' y'], '.-') 
[x y] = draw_dot(G<thres_conn);
gplot(G<thres_conn, [x' y'], '.-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 80;
V = 10*rand(n,10);%E_r=rand(n,2);E_r=E_r>(.4,.4);
E = ceil(n*rand(n,2));
[A,xy] = VxE_to_AdjacencyMatrix_and_XY_points(V,E);
spy(A);draw_dot(A);
gplot(A, xy, '.-'); 


n = 100;
p = 0.05; 
 G = rand(n,n) < p;
 G = triu(G,1);
 G = G + G';
thres_giant = 1/(n-1);
 [x y] = draw_dot(G);
gplot(G, [x' y'], '.-'); 



% CVX - GraphLaplacian Example 

% Absolute algebraic connectivity. Find edge weights that maximize the algebraic
% connectivity of the graph (i.e., the smallest positive eigenvalue of its Laplacian
% matrix). The optimal value is called the absolute algebraic connectivity by Fielder.
% 
% This generate a nice random graph - but no cut - n = 50; threshold = 0.2529;
% rand('state',209);
% xy = rand(n,2);
%A is an adjacency matrix
xy=zeros(N,2);
for i=1:N;
   for j=i:N;
     if A(i,j)>0;
      xy(i,1)=i;
      xy(i,2)=j;
     end;
   end;
 end;

% angle = 10*pi/180;
% Rotate = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
% xy = (Rotate*xy')';
% 
% Dist = zeros(n,n);
% for i=1:(n-1);
%   for j=i+1:n;
%     Dist(i,j) = norm( xy(i,:) - xy(j,:) );
%   end;
% end;
% Dist = Dist + Dist';
% Ad = Dist < threshold;
% Ad = Ad - eye(n);
% m = sum(sum(Ad))/2;
% 
% % find the incidence matrix
% A = zeros(n,m);
% l = 0;
% for i=1:(n-1);
%   for j=i+1:n;
%     if Ad(i,j)>0.5
%       l = l + 1;
%       A(i,l) =  1;
%       A(j,l) = -1;
%     end;
%   end;
% end;
% A = sparse(A);

% Compute edge weights: some optimal, some based on heuristics
[n,m] = size(A);


%Find xy




[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(A);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);

tau_fdla = 1/log(1/rho_fdla);
tau_fmmc = 1/log(1/rho_fmmc);
tau_md   = 1/log(1/rho_md);
tau_bc   = 1/log(1/rho_bc);
tau_mh   = 1/log(1/rho_mh);

% Plot results
figure(1), clf
plotgraph(A,xy,w_fdla);
text(0.55,1.05,'FDLA optimal weights')

figure(2), clf
plotgraph(A,xy,w_fmmc);
text(0.55,1.05,'FMMC optimal weights')

figure(3), clf
plotgraph(A,xy,w_md);
text(0.5,1.05,'Max degree optimal weights')

figure(4), clf
plotgraph(A,xy,w_bc);
text(0.5,1.05,'Best constant optimal weights')

figure(5), clf
plotgraph(A,xy,w_mh);
text(0.46,1.05,'Metropolis-Hastings optimal weights')


%$%%%%%%%%%%%%%%%%%%%%%%%% From Matlab BGL
%generate a second order finite difference approximation to the 
% Laplacian operator on a rectangular domain. This matrix does have a red-black ordering.
% n is the number of points used to discretize each dimension. 
% N is the total number rows and columns in the matrix/graph.
n = 12;
N = n*n;
A = delsq(numgrid('S',n+2));
%First, we compute a breadth first search on the graph and store the distance 
% each vertex is from the root. Because we really do not care, we'll choose 
% vertex 1 (row 1) of the matrix as the root vertex.
d = bfs(A,1);
d_even = find(mod(d,2) == 0);
d_odd = find(mod(d,2) == 1);%ptp =d_odd * d_even';
%To find the red-black ordering for an arbitrary matrix 
% (if we do not know it analytically) is easy using MatlabBGL.
% The key idea is to realize that a red-black ordering is 
% equivalent with the partition of vertices in a bipartite graph. 
% Once we see the problem in this light, we can quickly come up 
% with an algorithm that yields a potential red-black ordering.
% We begin by picking an arbitrary vertex and look at how far a 
% breadth first search goes at every step. To find the bipartition, 
% we look at all vertices which are an even distance from the root 
% and all the vertices which are an odd distance from the root. 
% If the matrix has a red-black ordering or is a bipartite graph, 
% this algorithm will find it.
% Implementing this algorithm is trivial using the MatlabBGL library. 
% First, we compute a breadth first search on the graph and store the 
% distance each vertex is from the root. Because we really do not care, 
% we'll choose vertex 1 (row 1) of the matrix as the root vertex. 
p = [d_odd' d_even'];
spy(A(p,p) - diag(diag(A(p,p))));
hold on;
plot([size(A,2)/2 size(A,2)/2],[0 size(A,1)], 'k-');
plot([0 size(A,2)],[size(A,1)/2 size(A,1)/2], 'k-');
hold off;
B=full(A);

C=zeros(N,N);
for i=1:(N);
  for j=i+1:N;
    if B(i,j)>=0
        C(i,j)=1;
        C(j,i)=1;
    end;
  end;
end;

draw_dot(C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOOST GRAPH LIB Examples

load ../graphs/max_flow_example.mat
max_flow(A,1,8)
[flow cut R F] = max_flow(A,1,8);
draw_dot(A);
full(R)     

load ../graphs/bfs_example.mat
[d dt pred] = bfs(A,2);
[ignore order] = sort(dt);
labels(order)
treeplot(pred);


load ../graphs/dfs_example.mat
[d dt ft pred] = dfs(A,2);
[ignore order] = sort(dt);
labels(order)
draw_dot(full(A));

    load graphs/clr-26-1.mat
    all_shortest_paths(A) ;B= full(A);
%  [D,P]
 D =  all_shortest_paths(A,struct('algname','johnson'));
 draw_dot(full(D));
%
    G1 = cycle_graph(9000,struct('directed',0));
    X1 = gursoy_atun_layout(G1,'topology','heart');
    G2 = grid_graph(50,50);
    X2 = gursoy_atun_layout(G2,'topology','square');
    G3 = grid_graph(50,50);
    X3 = gursoy_atun_layout(G3,'topology','circle');
    subplot(1,3,1); gplot(G1,X1,'k'); subplot(1,3,2); gplot(G2,X2,'k');
    subplot(1,3,3); gplot(G3,X3,'k');


As = cycle_graph(128,'directed',1); % compute a 6 node cycle graph
A = As; % set all the weights to be one initially
A(2,3) = 0; A(3,2) = 0; % make one edge have zero weight
fprintf('Edges\n');
full(As)
fprintf('Weights\n');
full(A)
X = kamada_kawai_spring_layout(G);
gplot(G,X,'.-');
% 
% Note that As is given as the graph in the following call, not A!
[d pred] = shortest_paths(As,1,'edge_weight',edge_weight_vector(As,A));
d(3) % distance from vertex 1 to vertex 3 should be just 1!

G = grid_graph(12,28);
X = kamada_kawai_spring_layout(G);
gplot(G,X,'.-');

G = grid_graph(6,5);
X = fruchterman_reingold_force_directed_layout(G);
gplot(G,X,'.-');

G = grid_graph(2*ones(1,8)); % compute 5d hypercube

G = grid_graph(2*ones(1,7)); % compute 5d hypercube
G = grid_graph(2*ones(1,6)); % compute 5d hypercube
G = grid_graph(2*ones(1,5)); % compute 5d hypercube
G = grid_graph(2*ones(1,4)); % compute 5d hypercube
G = grid_graph(2*ones(1,3)); % compute 5d hypercube
G = grid_graph(2*ones(1,2)); % compute 5d hypercube
G = grid_graph(2*ones(1,1)); % compute 5d hypercube


G = grid_graph(6,5);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');


% A grid in the xy plane is a planar graph.
G = grid_graph(6,5);
is_planar = boyer_myrvold_planarity_test(G)

% Recall that K_5 (the clique on 5 vertices) is not a planar graph.  Let's
% see what happens.
G = clique_graph(1);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(2);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(3);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(4);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(5);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(6);
X = gursoy_atun_layout(G);
gplot(G,X,'.-');

G = clique_graph(15);
X = fruchterman_reingold_force_directed_layout(G);
gplot(G,X,'.-');

is_planar = test_planar_graph(K5) % helpful wrapper

% We can also draw planar graphs
G = grid_graph(6,5);
X = chrobak_payne_straight_line_drawing(G);
gplot(G,X,'.-'); % it looks a little different!
% New option syntax
% You probably noticed that the "struct" command that permeated MatlabBGL
% calls before is gone in these examples.  We've moved to a new option
% syntax that gives you the _choice_ between the MatlabBGL struct style
% arguments and a list of key-value pairs

% We'll look at spanning trees on the clique graph with 5 vertices.  
% Using Prim's algorithm, the spanning tree we get depends on the root.  We
% always get a star graph rooted at the vertex we pick as the root.
G = clique_graph(5);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=j;x2=x1(:,:,1)/3+x1(:,:,2)/3+x1(:,:,3)/3;x3=double(x2)/255;
load earth % Load image data, X, and colormap, map
sphere; h = findobj('Type','surface');
hemisphere = [ones(789,931/2),... 
              x3,... 
              ones(789,931/2)];
set(h,'CData',flipud(hemisphere),'FaceColor','texturemap')
colormap(gray)
axis equal
view([90 0])
set(gca,'CameraViewAngleMode','manual')
view([65 30])

%Differential Calculus on a Matrix
n = 789;X=im2double(x1);
%Ml = load_image('lena',256);
X = load_image('1.jpg',256);X=x3;M=X;
colormap(gray);
M = rescale(crop(X,n));
options.bound = 'per';
T = compute_structure_tensor(M,1.5,8);
% Display the tensor fields. The color is proportional to the size of the tensor.
clf;
options.sub = 5;
plot_tensor_field(T, M, options);
% initial flow
v = perform_blurring(randn(n,n,2), 40, options);
[tmp,v] = compute_hodge_decompositon(v,options);
v = perform_vf_normalization(v);
% options for the PDE solver
dt = .3;
options.viscosity = 2*dt; % diffusion per frame
options.advspeed = 1*dt; % advection per frame
options.viscosity_texture = .3*dt; % diffusion of the texture
options.texture_histo = 'linear'; % fix the contrast
options.display = 0;
options.niter_fluid = 100;
% solve the PDE
[vlist,A] = perform_fluid_dynamics(v,M,options);
% display
sel = round( linspace(1,size(A,3),6) );
B = mat2cell(A(:,:,sel),n,n,ones(6,1));
clf;
imageplot(B{1});
imageplot(B{2});
imageplot(B{3});
imageplot(B{4});
imageplot(B{5});
imageplot(B{6});
imageplot(B{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three cluseters with a cut between each
clear all;
N=3*164;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.23;
p_G2 = 0.33;
p_G3 = 0.19;
p_Inbetween = 0.016;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;

[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));

B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween;
 A(G2, G3) = B_2 < p_Inbetween/3;
 A(G1, G3) = B_3 < p_Inbetween/5;

A = triu(A,1);
A = A + A';
spy(A);

L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);

spy(L);
spy(L2);
spy(L3-L2);

[V D] = eigs(L2, 3, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2)
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
plot(V(:,3), '.-');
plot(sort(V(:,3)), '.-');

spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

XY_ga = gursoy_atun_layout(sparse(A*1.0),'topology','circle');
XY_kk = kamada_kawai_spring_layout(sparse(A*1.0));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));
% wgplot(GEV_l,X);
% AB=sparse(double(Net));
% XB = gursoy_atun_layout(AB,'topology','heart');gplot(AB,X,'.-');

wgplot(A,XY_ga);
wgplot(A,XY_kk);
wgplot(A,XY_fr); %This is the best format for cut display
gplot(A,XY_fr);
WDEFAULT(gcf)
[row_clust_idx, col_clust_idx,y_index,x_index]=SpectralCoClustering(A,3);%,display,names)


[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Ic);
[ w_md   ] = max_deg(Ic);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);

Afmmc = abs(Ic*diag(w_fmmc)*Ic');
Amd = abs(Ic*diag(w_md)*Ic');
An = (An - diag(diag(An))) > 0.001;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));

%Generate sample from \beta Hermite ensemble
n=N;m=N;
beta=1;
d=sqrt(chi2rnd(beta*[n:-1:1]))';
H=spdiags(d,1,n,n)+spdiags(randn(n,1),0,n,n);
H=(H+H')/2;
[L Q]=eigs(H);
hist(L);v=diag(H);


gplot(An,XY_fr,'b-'); hold on;
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5)
gplot(A,XY_fr,'b:');
plot(xy(:,1), xy(:,2), 'ko','LineWidth',4, 'MarkerSize',4);
axis([0.05 1.1 -0.1 0.95]);
title('Subgraph with positive transition prob.')
hold off;


BB = triu(B_1,1);
BB = BB + BB';
 [W,S] = FMMC(BB) %gives a vector of the fastest mixing Markov chain
% edge weights for a graph described by the incidence matrix A (n x m).
% Here n is the number of nodes and m is the number of edges in the graph;
% each column of A has exactly one +1 and one -1.
%
% The FMMC edge weights are given the SDP:
%
%   minimize    s
%   subject to  -s*I <= I - L - (1/n)11' <= s*I
%               w >= 0,  diag(L) <= 1
%


%Maeshpart Pkg
specpart(A,XY_fr);
cutsize(A,ans);

geopart(A,XY_fr,2);
cutsize(A,ans);

metispart(sparse(A),sparse(XY_fr));
cutsize(A,ans)

gspart(A,XY_fr);
disp('cutsize(Tapir,ans)');
cutsize(A,ans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MESHPART

%%%%%%%%%%%%%%%%%%%%%%%%QUANTUM CLUSTERING - get the NN toolbox to run this
% % function [mapping] = qcClustering(sigma,steps,data)
% % performs quantum clustering + gradient descent
% % Input -
% % sigma - as in the algorithm
% % steps - number of steps
% % data - dims*samples matrix
% % doRecurse - if 1, then call the methods again - if not use the
% % bestCluster
% % numOfClusters - number of clusters
% % Output -
% % mapping results of the algorithm
%  [mapping] = qcClustering(2,20,XY_fr,1,1,3)









C_coeff=clustering_coefficients((sparse(A)));full(A)
hist(C_coeff,100);C_mu=mean(C_coeff);

Net = SFNG(N+20, 5, A);
C_coeff=clustering_coefficients(sparse(double(Net)));
hist(C_coeff,100);C_mu=mean(C_coeff);

[x y]=draw_dot(Net);
gplot(Net, [x' y'], '.-');












clear all;
N=150;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = randperm(N);
gs = N/2;
G1 = x(1:gs);
G2 = x(gs+1:end);
% Decide on the probabilities of edges within each group and between the two groups. 
% Because we are planting a partition, the probabilities of edges between the groups should be much 
% lower than the probability of edges within each group. Suppose that group 1 is a little more tightly connected than group 2.
p_G1 = 0.23;
p_G2 = 0.23;
p_Inbetween = 0.04;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(N-gs,N-gs) < p_G2;
A(G1, G2) = rand(gs, N-gs) < p_Inbetween;
At = triu(A,1);

A = At -At';
Ai=incidence(A);

Aa=incidence2adjacency(A);

L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);




[V D] = eigs(L2, 2, 'SA');
D(2,2);%ans = 46.7158
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:,2));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(At+At'*1.0));

wgplot(A,XY_fr); %This is the best format for cut display

xy=XY_fr;
%%%%%%%%%%%%%%%%%%%%%Meshpart Pkg
specpart(A,XY_fr);
cutsize(A,ans);

geopart(A,XY_fr,2);
cutsize(A,ans);

metispart(sparse(A),sparse(XY_fr));
cutsize(A,ans)

gspart(A,XY_fr);
disp('cutsize(Tapir,ans)');
cutsize(A,ans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MESHPART


[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Aa);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);


[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fdla);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fmmc);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_md);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_bc);
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_mh);

[he,hv]=wgPlot(A,XY_fr,'vertexWeight',w_fmmc,'vertexMetadata',w_fmmc*5,'vertexScale',200,'edgeColorMap',pink);
rand(length(A),1),'edgeColorMap',pink);











clear all;
memory;%clear, pack, whos, inmem, save, load, mlock, munlock
N=5*188;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/5;

G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3 = x(2*gs+1:3*gs);
G4 = x(3*gs+1:4*gs);
G5 = x(4*gs+1:5*gs);
%V= exp(1)./ log(1:5*N);plot(V);
V= (1)./ log(1:5*N);plot(V);
p_G1 = 0.53*V(1*N);
p_G2 = 0.43*V(1*N);
p_G3 = 0.39*V(1*N);
p_G4 = 0.29*V(1*N);
p_G5 = 0.19*V(1*N);

p_Inbetween_12 = 0.016*V(2*N);
p_Inbetween_13 = 0.010*V(2*N);
p_Inbetween_14 = 0.0108*V(2*N);
p_Inbetween_15 = 0.0106*V(2*N);

p_Inbetween_23 = 0.018*V(2*N);
p_Inbetween_24 = 0.0106*V(2*N);
p_Inbetween_25 = 0.0104*V(2*N);

p_Inbetween_34 = 0.018*V(2*N);
p_Inbetween_35 = 0.0106*V(2*N);

p_Inbetween_45 = 0.009*V(2*N);

A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;
A(G4, G4) = rand(gs,gs) < p_G4;
A(G5, G5) = rand(gs,gs) < p_G5;

B_1=rand(gs, gs);
B_2=rand(gs, gs);
B_3=rand(gs, gs);
B_4=rand(gs, gs);
B_5=rand(gs, gs);

A(G1, G2) = B_1 < p_Inbetween_12;
A(G1, G3) = B_1 < p_Inbetween_13;
A(G1, G4) = B_1 < p_Inbetween_14;
A(G1, G5) = B_1 < p_Inbetween_15;


A(G2, G3) = B_2 < p_Inbetween_23;
A(G2, G4) = B_2 < p_Inbetween_24;
A(G2, G5) = B_2 < p_Inbetween_25;


A(G3, G4) = B_3 < p_Inbetween_34;
A(G3, G5) = B_3 < p_Inbetween_35;

A(G4, G5) = B_4 < p_Inbetween_45;

A = triu(A,1);
A = A + A';
spy(A);

%L = del2(A);
%L2=laplacian(A);%

%Metis Test
nparts=5;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));
XY = kamada_kawai_spring_layout(sparse(A*1.0));

gplot(A,XY_fr);
map = metisdice(sparse(A),nparts,XY_fr);


Ic = adjacency2incidence(A);
L3 = Ic*Ic';%compute_combinatorial_laplacian(A);
[V D] = eigs(L3, 5, 'SA');
sigma=eig(A);
hist(sigma,50);
D(2,2)

plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');

plot(V(:,3), '.-');
plot(sort(V(:,3)), '.-');

plot(V(:,4), '.-');
plot(sort(V(:,4)), '.-');

plot(V(:,5), '.-');
plot(sort(V(:,5)), '.-');


spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));


wgplot(A,XY_fr); %This is the best format for cut display
gplot(A,XY_fr);

memory

C_coeff=clustering_coefficients(sparse(A));
plot(hist(C_coeff,100));

[ w_fdla ] = fdla(sparse(A));
[ w_fmmc ] = fmmc(Ic);
[ w_md   ] = max_deg(Ic);
[ w_bc   ] = best_const(A);
[ w_mh   ] = mh(A);

Afmmc = abs(Ic*diag(w_fmmc)*Ic');
%Amd = abs(Ic*diag(w_md)*Ic');
% Comment this out to keep weighted adjacency matrix

An = (An - diag(diag(An))) > 0.001;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));


%Waxman random network topology generator
%The random network topology generator introduced in Waxman (1988) is a geographic model for the growth of a computer network. In this model the nodes of the network are uniformly distributed in the plane and edges are added according to probabilities that depend on the distances between the nodes. The probability to have an edge between nodes u and v is given by 
%P(u, v) = ae-d/(bL)  
lambda=2,alpha=1, beta=.5, domain=[-1 1 -1 1] .* 3;
[adj_matr, nd_coord] = waxtop(lambda,alpha, beta, domain); 
C_coeff=clustering_coefficients(sparse(A));
plot(hist(C_coeff,100));

% Inputs: 
%   lambda - intensity of the Poisson process
%   alpha - maximal link probability
%   beta - parameter to control length of the edges. Increased <beta>
%     yields a larger ratio of long edges to short edges
%   domain - bounds for the region. A 4-dimensional vector in
%     the form [x_min x_max y_min y_max]


% %Galton-Watson 
ngen=7;
p=[0 1 1 1 2 3 4 3 4 3 4 6 7 8 10 10 11 11 11 12 15 16 18 17 17 17 21 23 24 24 22 25 22 25 22 25 27 28 34 31 31 29 30 36 29 30 36 29 30 36];
[parents] = branch(ngen, p);

% There was concern amongst the Victorians that aristocratic surnames were becoming extinct. Galton originally posed the question regarding the probability of such an event in the Educational Times of 1873, and the Reverend Henry William Watson replied with a solution. Together, they then wrote an 1874 paper entitled On the probability of extinction of families. Galton and Watson appear to have derived their process independently of the earlier work by I. J. Bienaymé; see Heyde and Seneta 1977. For a detailed history see Kendall (1966 and 1975).
% 
% 
% [edit] Concepts
% Assume, as was taken for granted in Galton's time, that surnames are passed on to all male children by their father. Suppose the number of a man's sons to be a random variable distributed on the set { 0, 1, 2, 3, ...}. Further suppose the numbers of different men's sons to be independent random variables, all having the same distribution.
% 
% Then the simplest substantial mathematical conclusion is that if the average number of a man's sons is 1 or less, then their surname will surely die out, and if it is more than 1, then there is more than zero probability that it will survive forever.
% 
% Modern applications include the survival probabilities for a new mutant gene, or the initiation of a nuclear chain reaction, or the dynamics of disease outbreaks in their first generations of spread, or the chances of extinction of small population of organisms; as well as explaining (perhaps closest to Galton's original interest) why only a handful of males in the deep past of humanity now have any surviving male-line descendants, reflected in a rather small number of distinctive human Y-chromosome DNA haplogroups.
% 
% A corollary of high extinction probabilities is that if a lineage has survived, it is likely to have experienced, purely by chance, an unusually high growth rate in its early generations at least when compared to the rest of the population.
% 
% 
% [edit] Mathematical definition
% A Galton-Watson process is a stochastic process {Xn} which evolves according to the recurrence formula X0 = 1 and
% 
%  
% where for each n,  is a sequence of IID natural number-valued random variables. The extinction probability (i.e. the probability of final extinction)is given by
% 
%  
% This is clearly equal to zero if each member of the population has exactly one descendent. Excluding this case (usually called the trivial case) there exists a simple necessary and sufficient condition, which is:
% 
% 
% [edit] Extinction criterion for Galton-Watson process
% In the non-trivial case the probability of final extinction is equal to one if E{?1} ? 1 and strictly less than one if E{?1} > 1.
% 
% The process can be treated analytically using the method of probability generating functions.
% 
% If the number of children ? j at each node follows a Poisson distribution, a particularly simple recurrence can be found for the total extinction probability xn for a process starting with a single individual at time n = 0:
% 
%  
% giving the curves plotted above.
% 
% 
% [edit] Bisexual Galton–Watson process
% In the (classical) Galton–Watson process defined above, only men count, that is, the reproduction can be understood as being asexual. The more natural corresponding version for (bi)sexual reproduction is the so-called 'bisexual Galton–Watson process', where only couples can reproduce. In this process, each child is supposed to be male or female, independently of each other, with a specified probability, and a so-called 'mating function' determines how many couples will form in a given generation. As before, reproduction of different couples are considered to be independent of each other. Now the analogue of the trivial case corresponds to the case of each male and female reproducing in exactly one couple, having one male and one female descendent, and that the mating function takes the value of the minimum of the number of males and females (which are then the same from the next generation onwards).
% 
% Since the total reproduction within a generation depends now strongly on the mating function, there exists in general no simple necessary and sufficient condition for final extinction as it is the case in the classical Galton–Watson process. However, excluding the non-trivial case, the concept of the averaged reproduction mean (Bruss (1984)) allows for a general sufficient condition for final extinction, treated in the next section.
% 
% 
% [edit] Extinction criterion (bisexual Galton-Watson process)
% If in the non-trivial case the averaged reproduction mean per couple stays bounded over all generations and will not exceed 1 for a sufficiently large population size, then the probability of final extinction is always 1.
% 
% 
% [edit] Examples
% Countries that have used family names for many generations exhibit the Galton–Watson process in their low number of surviving family names:
% 
% Korean names are the most striking example, with 250 family names, and 45% of the population sharing three family names 
% Chinese names are similar, with 22% of the population sharing three family names (numbering close to 300 million people), and the top 200 names covering 96% of the population. 
% By contrast:
% 
% Dutch names have only included a family name since the Napoleonic Wars in the early 19th century, and there are over 68,000 Dutch family names. 
% Thai names have only included a family name since 1920, and only a single family can use a given family name, hence there are a great number of Thai names. Further, Thai people change their family names with some frequency, complicating the analysis. 
% 
% [edit] See also
% Branching process 
% 
% [edit] References
% F. Thomas Bruss (1984). "A Note on Extinction Criteria for Bisexual Galton–Watson Processes". Journal of Applied Probability 21: 915–919. 
% C C Heyde and E Seneta (1977). I.J. Bienayme: Statistical Theory Anticipated. Berlin, Germany. 
% D G Kendall (1966). Journal of the London Mathematical Society 41: 385–406 
% D G Kendall (1975). Bulletin of the London Mathematical Society 7: 225–253 
% H W Watson and Francis Galton, "On the Probability of the Extinction of Families", Journal of the Anthropological Institute of Great Britain, volume 4, pages 138–144, 1875. 
% 
% [edit] External links
% The original Galton–Watson paper: On the Probability of the Extinction of Families 
% "Survival of a Single Mutant" by Peter M. Lee of the University of York 
% Retrieved from "http://en.wikipedia.org/wiki/Galton%E2%80%93Watson_process"

%%%%%%%%%%%%%%%%%%%%%%RMT
n=64;
X=rand(n);
A=(X+X')/2;
[Q,L]=eig(A);
L=diag(L);
J=zeros(n * (n+1)/2);
epsilon=1e-7;
idx=1
mask=triu(ones(n),1);
mask=logical(mask(:,:));
for i=1:n
    for j=1:n
        E_ij=zeros(n);
        E_ij(i,j)=1;
        E_ij(j,i)=1;
        Ap=A+epsilon*E_ij;
        [Qp,Lp]=eig(Ap);
        dL=(diag(Lp)-L)/epsilon;%eigenvalue peturbation
        QdQ=Q'*(Qp-Q)*epsilon;%eigenvector peturbation
        J(1:n,idx)=dL;% Eigenvalue portion 
        J(n+1:end,idx)=QdQ(mask);
        idx=idx+1;
    end
end

diff=(1/abs(det(J))-1/abs(vander(L)));

%Generate sample from \beta Hermite ensemble

m=n;
beta=1;
d=sqrt(chi2rnd(beta*[n:-1:1]))';
H=spdiags(d,1,n,n)+spdiags(randn(n,1),0,n,n);
H=(H+H')/2;
[L Q]=eigs(H);
hist(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Call some functions in RMT
n=1024;m=1024;isreal=1;
W=wishart(n,m,isreal);
Wi=wigner(n,isreal);
t = trijacobi(n,a,b);
maxmom=20;
m = momwigner(n,maxmom,isreal);
f = manovaedf(c1,c2);
m1=512;m2=512;
j = manova(n,m1,m2,isreal);C=eig(j)
plot(hist(C,200));
js = jsymeig(a);%JSYMEIG  Jacobian for the symmetric eigenvalue problem
[h,hn,xspan]=histn(data,x0,binsize,xf);%HISTN Normalized Histogram.
%    [H,HN,XSPAN] = HISTN(DATA,X0,BINSIZE,XF) generates the normalized
%    histogram of area 1 from the values in DATA which are binned into
%    equally spaced containers that span the region from X0 to XF
%    with a bin width specified by BINSIZE.

N=1024%4096;
%Generate a Winger - GOE
G = randn(N)/sqrt(N); 
A = (G+G')/sqrt(2); 
C1=eig(A);
C_Epdf1=hist(C1,200);
plot(C_Epdf1);

C1=eig(W);
C_Epdf1=hist(C1,200);
plot(C_Epdf1);


%Generate Wishart Matrix
G = randn(N,2*N)/sqrt(2*N); 
B = G*G';
C=eig(B);
EmpiricalPDF= hist(C,200);
plot(EmpiricalPDF);
N=2048%4096;
%Generate a Winger - GOE
G = randn(N)/sqrt(N); 
A = (G+G')/sqrt(2); 
C1=eig(A);
C_Epdf1=hist(C1,200);

