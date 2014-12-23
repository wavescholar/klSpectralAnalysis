%051709
clear all;
N=3*    44;
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

A = A + A';A=sparse(A);spy(A)
xy=fruchterman_reingold_force_directed_layout(max(A,A'),'force_pairs','grid','initial_temp',1000,'iteration',1); 
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));


A(find(A))=(rand(size(find(A)))+0.01)*100;  % make non-zero element of A random.
figure;
tic;
[he,hv]=wgPlot(A,XY_fr,'vertexWeight',(1:N),'vertexMetadata',(1:N),'edgeColorMap',pink);
toc

[x y]=draw_dot(A);

[hE,hV]=wgPlot(A,coord,varargin)


A(G1, G1) = A(G1, G1) < p_G1;
A(G2, G2) = A(G2, G2) < p_G2;
A(G3, G3) = A(G3, G3) < p_G3;
B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween;
 A(G2, G3) = B_2 < p_Inbetween/3;
 A(G1, G3) = B_3 < p_Inbetween/5;
 A = A + A';
 
 
 
[row_clust_idx, col_clust_idx,y_index,x_index]= SpectralCoClustering(A,3);





 B = allspath(A);










[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));



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

