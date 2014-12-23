N=4096;

%Generate a Winger - GOE
G = randn(N)/sqrt(N); 
A = (G+G')/sqrt(2); %Winger Matrix
imagesc(A);I=round(A*255+200);
imagesc(I,[0 255]); colormap(gray);
C1=eig(A);
hist(C1,100);
C_Epdf1=hist(C1,200);
%Sample again from GOE
G = randn(N)/sqrt(N); 
A = (G+G')/sqrt(2);
C2=eig(A);
hist(C2,100);
C_Epdf2=hist(C2,200);
%Calculate KL Divergence of E_pdfs
kl = kldivergence (C1,C2);

g=C_Epdf1-C_Epdf2;
%Generate Wishart Matrix
G = randn(N,2*N)/sqrt(2*N); 
B = G*G';
%RMtool from rao requires the symbolic toolbox
% syms m z
% LmzA = numden(m-0.5/(1-z)-0.5/(2-z));
% LmzB = AtimesWish(LmzA,0.5);
C=eig(B);
max_minus_min_eig=range(C);
EmpiricalPDF= hist(B,200);
plot(EmpiricalPDF);
I=round(B*255+200);
imagesc(I,[0 255]); colormap(gray);

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

At = At -At';

[ w_fdla ] = fdla(sparse(At));
[ w_fmmc ] = fmmc(Net);
[ w_md   ] = max_deg(A);
[ w_bc   ] = best_const(Net);
[ w_mh   ] = mh(Net);


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
L = del2(A);%L=compute_laplacian(A);%L=compute_combinatorial_laplacian(A);
[V D] = eigs(L, 2, 'SA');
D(2,2);%ans = 46.7158
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');
[ignore p] = sort(V(:,2));
spy(A(p,p));draw_dot(A(p,p));
%[ignore p] = sort(V(:,2));
%spy(A(p,p));
%Let's do an MDS on the graph adjacency
 [points,vaf]=mds(A,4);%metric,iterations,learnrate)
 %D=all_shortest_paths(A);
 plot(points(:,1),points(:,2));
[ w_mh   ] = mh(A);
plotgraphBruce(A,points,w_mh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Generation of cluster grapg [my own term ]
% Three cluseters with a cut between each
clear all;
N=3*128;
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
 A(G3, G1) = B_3 < p_Inbetween/50;

A = triu(A,1);
A = A + A';spy(A);
L = del2(A);spy(L);
[V D] = eigs(L, 2, 'SA');sigma=eig(A);hist(sigma,50);
D(2,2)
plot(V(:,2), '.-');
plot(sort(V(:,2)), '.-');spy(A);%draw_dot(A);
[ignore p] = sort(V(:,2));
spy(A(p,p));

[x y]=draw_dot(double(A));
gplot(A, [x' y'], '.-')

C_coeff=clustering_coefficients((sparse(A)));full(A)
hist(C_coeff,100);C_mu=mean(C_coeff);

Net = SFNG(N+20, 5, A);
C_coeff=clustering_coefficients(sparse(double(Net)));
hist(C_coeff,100);C_mu=mean(C_coeff);

[x y]=draw_dot(Net);
gplot(Net, [x' y'], '.-');


    sA=sparse(A);
    X1 = gursoy_atun_layout(sA,'topology','heart');gplot(sA,X1,'k');
    

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
[cn csz]=core_numbers(A);  cs=unique(cn);spy(A);
 %[A,xy] = VxE_to_AdjacencyMatrix_and_XY_points(V,E);
 X = gursoy_atun_layout(A,'topology','heart');
 %X = kamada_kawai_spring_layout(A);%fruchterman_reingold_force_directed_layout(A);
gplot(A,X,'.-');


C_coeff=clustering_coefficients(A)
hist(C_coeff,100);

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
    G1 = cycle_graph(5000,struct('directed',0));
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


%Differential Calculus on a Matrix
n = 256;
%Ml = load_image('lena',256);
M = load_image('0409061117g',256);
%colormap(gray);
M = rescale(crop(M,n));
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
imageplot(B);


