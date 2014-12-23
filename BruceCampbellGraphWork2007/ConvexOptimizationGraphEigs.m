clear all;
N=3*64;
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
% spy(A);
% 
% L = del2(A);
L2=laplacian(A);%

Ic = adjacency2incidence(A);
% L3 = Ic*Ic';%compute_combinatorial_laplacian(A);

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

% XY_ga = gursoy_atun_layout(sparse(A*1.0),'topology','heart');
% XY_kk = kamada_kawai_spring_layout(sparse(A*1.0));
% XY_fr=fruchterman_reingold_force_directed_layout(sparse(A*1.0));
% % wgplot(GEV_l,X);
% % AB=sparse(double(Net));
% % XB = gursoy_atun_layout(AB,'topology','heart');gplot(AB,X,'.-');
% 
% wgplot(A,XY_ga);
% wgplot(A,XY_kk);
% wgplot(A,XY_fr); %This is the best format for cut display
% gplot(A,XY_fr);
% WDEFAULT(gcf)
% [row_clust_idx, col_clust_idx,y_index,x_index]=SpectralCoClustering(A,3);%,display,names)


[ w_fdla ] = fdla(Ic);
[ w_fmmc ] = fmmc(Ic);
[ w_md   ] = max_deg(Ic);
[ w_bc   ] = best_const(Ic);
[ w_mh   ] = mh(Ic);

Afmmc = abs(Ic*diag(w_fmmc)*Ic');
Amd = abs(Ic*diag(w_md)*Ic');
An = (An - diag(diag(An))) > 0.001;
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Afmmc*1.0));
wgplot(Afmmc,XY_fr,'vertexWeight',V(:,3));

XY_kk=kamada_kawai_spring_layout(sparse(Afmmc));
XY_ga=gursoy_atun_layout(sparse(Afmmc));
XY_fr=fruchterman_reingold_force_directed_layout(sparse(Amd*1.0));
wgplot(Amd,XY_fr,'vertexWeight',V(:,3));


XY_fr=fruchterman_reingold_force_directed_layout(sparse(A));
wgplot(A,XY_kk,'vertexWeight',w_fdla);

cm=[w_fmmc,w_fmmc,w_fmmc];
tic;
[he,hv]=wgPlot(Afmmc,XY_kk,'vertexWeight',V(:,3),'edgecolormap',pink);
toc