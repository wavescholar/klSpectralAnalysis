% Three cluseters with a cut between each
clear all;
N=3*228;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/3;
G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3= x(2*gs+1:3*gs);
p_G1 = 0.02;
p_G2 = 0.05;
p_G3 = 0.03;
p_Inbetween = 0.004;
A(G1, G1) = rand(gs,gs) < p_G1;
A(G2, G2) = rand(gs,gs) < p_G2;
A(G3, G3) = rand(gs,gs) < p_G3;

[n1,m1]=size( A(G1, G2));
[n2,m2]=size( A(G2, G3));
[n3,m3]=size( A(G3, G1));

B_1=rand(gs, gs);B_2=rand(gs, gs);B_3=rand(gs, gs);
 A(G1, G2) = B_1 < p_Inbetween/7;
 A(G2, G3) = B_2 < p_Inbetween/5;
 A(G1, G3) = B_3 < p_Inbetween/9;

A = triu(A,1);
A = A + A';
spy(A);  

%Metis Test
nparts=3;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0),'initial_temp',1000);
wgplot(A,XY_fr);
map = metisdice(sparse(A),nparts,XY_fr);
K=sparse(A);save('dunno.txt','K', '-double','-ASCII');
csvwrite('dunno4.txt',A)
fid = fopen('dunno3.txt', 'wt');
[n,m]=size(A);
for i =1:m
fprintf(fid, '%d\n', A(i,:));
end 
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
N=5*32;
% In this segment, we'll plant a partition in a graph, and then use the second smallest eigenvector to find it.
% As always, the first step is to generate our dataset. In this example, we'll be a little more ambitious and use a larger number of vertices. 
x = (1:N);%randperm(N);
gs = N/5;

G1 = x(1:gs);
G2 = x(gs+1:2*gs);
G3 = x(2*gs+1:3*gs);
G4 = x(3*gs+1:4*gs);
G5 = x(4*gs+1:5*gs);

p_G1 = 0.53;
p_G2 = 0.43;
p_G3 = 0.39;
p_G4 = 0.29;
p_G5 = 0.19;

p_Inbetween_12 = 0.016;
p_Inbetween_13 = 0.010;
p_Inbetween_14 = 0.0108;
p_Inbetween_15 = 0.0106;

p_Inbetween_23 = 0.018;
p_Inbetween_24 = 0.0106;
p_Inbetween_25 = 0.0104;

p_Inbetween_34 = 0.018;
p_Inbetween_35 = 0.0106;

p_Inbetween_45 = 0.016;

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